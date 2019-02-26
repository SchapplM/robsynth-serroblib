% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:30
% EndTime: 2019-02-26 20:25:30
% DurationCPUTime: 0.12s
% Computational Cost: add. (178->30), mult. (188->43), div. (0->0), fcn. (141->9), ass. (0->19)
t183 = pkin(10) + qJ(4);
t179 = sin(t183);
t181 = cos(t183);
t185 = sin(pkin(11));
t186 = cos(pkin(11));
t193 = r_i_i_C(1) * t186 - r_i_i_C(2) * t185 + pkin(4);
t198 = r_i_i_C(3) + qJ(5);
t199 = t193 * t179 - t198 * t181;
t200 = t199 * qJD(4) - t179 * qJD(5);
t184 = qJ(1) + pkin(9);
t180 = sin(t184);
t197 = qJD(1) * t180;
t182 = cos(t184);
t196 = qJD(1) * t182;
t195 = qJD(4) * t182;
t192 = t185 * r_i_i_C(1) + t186 * r_i_i_C(2) + pkin(7) + qJ(3);
t191 = -t198 * t179 - t193 * t181;
t189 = -cos(pkin(10)) * pkin(3) - pkin(2) + t191;
t1 = [t182 * qJD(3) + t200 * t180 + (-cos(qJ(1)) * pkin(1) - t192 * t180 + t189 * t182) * qJD(1), 0, t196 (t193 * t197 - t198 * t195) * t179 + (-t198 * t197 + (-t193 * qJD(4) + qJD(5)) * t182) * t181, -t179 * t197 + t181 * t195, 0; t180 * qJD(3) - t200 * t182 + (-sin(qJ(1)) * pkin(1) + t192 * t182 + t189 * t180) * qJD(1), 0, t197, -t199 * t196 + (t191 * qJD(4) + qJD(5) * t181) * t180, t180 * qJD(4) * t181 + t179 * t196, 0; 0, 0, 0, -t200, qJD(4) * t179, 0;];
JaD_transl  = t1;
