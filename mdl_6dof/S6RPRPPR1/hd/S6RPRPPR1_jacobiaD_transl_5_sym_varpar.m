% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:12
% EndTime: 2019-02-26 20:39:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (191->28), mult. (208->38), div. (0->0), fcn. (154->10), ass. (0->20)
t186 = qJ(3) + pkin(10);
t182 = sin(t186);
t184 = cos(t186);
t188 = sin(pkin(11));
t189 = cos(pkin(11));
t200 = r_i_i_C(1) * t189 - r_i_i_C(2) * t188 + pkin(4);
t205 = r_i_i_C(3) + qJ(5);
t195 = t200 * t182 - t205 * t184 + sin(qJ(3)) * pkin(3);
t209 = t195 * qJD(3) - t182 * qJD(5);
t208 = -t205 * t182 - t200 * t184 - cos(qJ(3)) * pkin(3);
t187 = qJ(1) + pkin(9);
t183 = sin(t187);
t204 = qJD(1) * t183;
t185 = cos(t187);
t203 = qJD(1) * t185;
t202 = qJD(3) * t184;
t199 = t188 * r_i_i_C(1) + t189 * r_i_i_C(2) + pkin(7) + qJ(4);
t197 = -pkin(2) + t208;
t194 = t208 * qJD(3) + qJD(5) * t184;
t1 = [t185 * qJD(4) + t209 * t183 + (-cos(qJ(1)) * pkin(1) - t199 * t183 + t197 * t185) * qJD(1), 0, t194 * t185 + t195 * t204, t203, -t182 * t204 + t185 * t202, 0; t183 * qJD(4) - t209 * t185 + (-sin(qJ(1)) * pkin(1) + t199 * t185 + t197 * t183) * qJD(1), 0, t194 * t183 - t195 * t203, t204, t182 * t203 + t183 * t202, 0; 0, 0, -t209, 0, qJD(3) * t182, 0;];
JaD_transl  = t1;
