% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:45
% EndTime: 2019-02-26 20:39:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (161->25), mult. (164->34), div. (0->0), fcn. (113->8), ass. (0->18)
t159 = qJ(3) + pkin(10);
t155 = sin(t159);
t157 = cos(t159);
t174 = r_i_i_C(3) + qJ(5);
t178 = pkin(4) - r_i_i_C(2);
t166 = t178 * t155 - t174 * t157 + sin(qJ(3)) * pkin(3);
t180 = t166 * qJD(3) - t155 * qJD(5);
t179 = -t174 * t155 - t178 * t157 - cos(qJ(3)) * pkin(3);
t175 = r_i_i_C(1) + qJ(4) + pkin(7);
t160 = qJ(1) + pkin(9);
t156 = sin(t160);
t173 = qJD(1) * t156;
t158 = cos(t160);
t172 = qJD(1) * t158;
t171 = qJD(3) * t157;
t168 = -pkin(2) + t179;
t165 = t179 * qJD(3) + qJD(5) * t157;
t1 = [t158 * qJD(4) + t180 * t156 + (-cos(qJ(1)) * pkin(1) - t175 * t156 + t168 * t158) * qJD(1), 0, t165 * t158 + t166 * t173, t172, -t155 * t173 + t158 * t171, 0; t156 * qJD(4) - t180 * t158 + (-sin(qJ(1)) * pkin(1) + t175 * t158 + t168 * t156) * qJD(1), 0, t165 * t156 - t166 * t172, t173, t155 * t172 + t156 * t171, 0; 0, 0, -t180, 0, qJD(3) * t155, 0;];
JaD_transl  = t1;
