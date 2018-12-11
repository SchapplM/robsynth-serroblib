% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_2_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_2_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_2_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_2_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_2_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.13s
% Computational Cost: add. (89->31), mult. (138->53), div. (0->0), fcn. (94->9), ass. (0->27)
t192 = sin(pkin(6)) * (pkin(10) + r_i_i_C(3));
t173 = pkin(6) + qJ(2);
t190 = cos(t173) / 0.2e1;
t174 = pkin(6) - qJ(2);
t189 = cos(t174);
t188 = sin(t173);
t187 = -qJD(2) / 0.2e1;
t171 = sin(t174);
t186 = qJD(2) * t171;
t177 = sin(qJ(1));
t185 = qJD(2) * t177;
t179 = cos(qJ(1));
t184 = qJD(2) * t179;
t164 = t189 / 0.2e1 + t190;
t162 = t164 * qJD(2);
t176 = sin(qJ(2));
t183 = -t162 * t179 + t176 * t185;
t163 = t188 / 0.2e1 - t171 / 0.2e1;
t178 = cos(qJ(2));
t182 = t163 * t179 + t177 * t178;
t181 = t163 * t177 - t178 * t179;
t180 = -t177 * t162 - t176 * t184;
t165 = t188 * t187;
t161 = t186 / 0.2e1 + t165;
t160 = t178 * t185 - t161 * t179 + (t164 * t177 + t176 * t179) * qJD(1);
t159 = -t178 * t184 - t177 * t161 + (-t164 * t179 + t176 * t177) * qJD(1);
t1 = [t183 * r_i_i_C(1) + t160 * r_i_i_C(2) + (-t179 * pkin(1) + t181 * r_i_i_C(1) - t177 * t192) * qJD(1), t159 * r_i_i_C(1) + (t182 * qJD(1) - t180) * r_i_i_C(2), 0, 0, 0, 0; t180 * r_i_i_C(1) + t159 * r_i_i_C(2) + (-t177 * pkin(1) - t182 * r_i_i_C(1) + t179 * t192) * qJD(1), -t160 * r_i_i_C(1) + (t181 * qJD(1) + t183) * r_i_i_C(2), 0, 0, 0, 0; 0 (qJD(2) * t190 + t189 * t187) * r_i_i_C(1) + (t165 - t186 / 0.2e1) * r_i_i_C(2), 0, 0, 0, 0;];
JaD_transl  = t1;
