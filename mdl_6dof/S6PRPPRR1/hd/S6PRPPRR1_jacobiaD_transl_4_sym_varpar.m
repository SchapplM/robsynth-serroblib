% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->24), mult. (241->43), div. (0->0), fcn. (228->10), ass. (0->23)
t203 = -r_i_i_C(3) - qJ(4);
t202 = qJD(2) * pkin(2);
t196 = cos(pkin(6));
t197 = sin(qJ(2));
t201 = t196 * t197;
t190 = sin(pkin(11));
t194 = cos(pkin(11));
t198 = cos(qJ(2));
t200 = t190 * t198 + t194 * t197;
t188 = t197 * t190 - t198 * t194;
t199 = cos(pkin(12)) * r_i_i_C(1) - sin(pkin(12)) * r_i_i_C(2) + pkin(3);
t185 = t200 * t196;
t187 = t200 * qJD(2);
t186 = t188 * qJD(2);
t195 = cos(pkin(10));
t192 = sin(pkin(6));
t191 = sin(pkin(10));
t184 = qJD(2) * t185;
t183 = t196 * t186;
t182 = t192 * t187;
t179 = t191 * t184 + t195 * t186;
t177 = -t195 * t184 + t191 * t186;
t1 = [0 -(t191 * t185 + t195 * t188) * qJD(4) + t203 * (-t191 * t183 + t195 * t187) + (t191 * t201 - t195 * t198) * t202 + t199 * t179, 0, -t179, 0, 0; 0 -(-t195 * t185 + t191 * t188) * qJD(4) + t203 * (t195 * t183 + t191 * t187) + (-t191 * t198 - t195 * t201) * t202 + t199 * t177, 0, -t177, 0, 0; 0, -t199 * t182 + (t200 * qJD(4) + t203 * t186 - t197 * t202) * t192, 0, t182, 0, 0;];
JaD_transl  = t1;
