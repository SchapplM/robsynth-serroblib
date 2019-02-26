% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:59
% EndTime: 2019-02-26 20:41:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (128->25), mult. (216->33), div. (0->0), fcn. (160->8), ass. (0->20)
t174 = qJ(3) + pkin(9);
t171 = sin(t174);
t172 = cos(t174);
t175 = sin(pkin(10));
t176 = cos(pkin(10));
t187 = r_i_i_C(1) * t176 - r_i_i_C(2) * t175 + pkin(4);
t192 = r_i_i_C(3) + qJ(5);
t185 = t187 * t171 - t192 * t172 + sin(qJ(3)) * pkin(3);
t199 = qJD(4) + (qJ(2) + t185) * qJD(1);
t197 = t185 * qJD(3) - qJD(5) * t171;
t184 = t192 * t171 + t187 * t172 + cos(qJ(3)) * pkin(3);
t196 = -t184 * qJD(3) + t172 * qJD(5);
t179 = sin(qJ(1));
t191 = qJD(1) * t179;
t181 = cos(qJ(1));
t173 = qJD(1) * t181;
t190 = qJD(3) * t171;
t183 = qJD(1) * t184;
t182 = qJD(2) + (-t175 * r_i_i_C(1) - t176 * r_i_i_C(2) - pkin(1) - pkin(7) - qJ(4)) * qJD(1) - t196;
t1 = [-t199 * t179 + t182 * t181, t173, -t197 * t179 + t181 * t183, -t191, -t172 * t173 + t179 * t190, 0; t182 * t179 + t199 * t181, t191, t179 * t183 + t197 * t181, t173, -t172 * t191 - t181 * t190, 0; 0, 0, t196, 0, qJD(3) * t172, 0;];
JaD_transl  = t1;
