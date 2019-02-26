% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:25
% EndTime: 2019-02-26 20:41:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (153->26), mult. (232->39), div. (0->0), fcn. (176->7), ass. (0->20)
t171 = sin(pkin(10));
t172 = cos(pkin(10));
t170 = pkin(9) + qJ(3);
t168 = sin(t170);
t169 = cos(t170);
t183 = r_i_i_C(1) * t171 + r_i_i_C(2) * t172 + qJ(4);
t185 = pkin(3) + r_i_i_C(3) + qJ(5);
t180 = t185 * t168 - t183 * t169;
t176 = -t180 * qJD(3) + t168 * qJD(4) + t169 * qJD(5);
t191 = t176 + (t172 * r_i_i_C(1) - t171 * r_i_i_C(2) + pkin(4) + pkin(7) + qJ(2)) * qJD(1);
t174 = sin(qJ(1));
t189 = qJD(1) * t174;
t175 = cos(qJ(1));
t188 = qJD(1) * t175;
t187 = qJD(3) * t174;
t186 = qJD(3) * t175;
t181 = -t183 * t168 - t185 * t169;
t178 = qJD(2) + (-cos(pkin(9)) * pkin(2) - pkin(1) + t181) * qJD(1);
t177 = t181 * qJD(3) + qJD(4) * t169 - qJD(5) * t168;
t1 = [-t191 * t174 + t178 * t175, t188, t177 * t175 + t180 * t189, -t168 * t189 + t169 * t186, -t168 * t186 - t169 * t189, 0; t178 * t174 + t191 * t175, t189, t177 * t174 - t180 * t188, t168 * t188 + t169 * t187, -t168 * t187 + t169 * t188, 0; 0, 0, t176, qJD(3) * t168, qJD(3) * t169, 0;];
JaD_transl  = t1;
