% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:02
	% EndTime: 2019-12-29 18:17:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:07
	% EndTime: 2019-12-29 18:17:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:02
	% EndTime: 2019-12-29 18:17:02
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:09
	% EndTime: 2019-12-29 18:17:09
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(8);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(6);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0; 0, -t33, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:09
	% EndTime: 2019-12-29 18:17:09
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->20), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->17)
	t145 = qJ(2) + pkin(8);
	t143 = sin(t145);
	t144 = cos(t145);
	t161 = r_i_i_C(3) + qJ(4);
	t165 = pkin(3) + r_i_i_C(1);
	t154 = t165 * t143 - t161 * t144 + sin(qJ(2)) * pkin(2);
	t151 = -t154 * qJD(2) + t143 * qJD(4);
	t168 = t151 + (r_i_i_C(2) + qJ(3) + pkin(6)) * qJD(1);
	t167 = -t161 * t143 - t165 * t144 - cos(qJ(2)) * pkin(2);
	t148 = sin(qJ(1));
	t160 = qJD(1) * t148;
	t150 = cos(qJ(1));
	t159 = qJD(1) * t150;
	t158 = qJD(2) * t144;
	t153 = qJD(3) + (-pkin(1) + t167) * qJD(1);
	t152 = t167 * qJD(2) + qJD(4) * t144;
	t1 = [-t168 * t148 + t153 * t150, t152 * t150 + t154 * t160, t159, -t143 * t160 + t150 * t158, 0; t153 * t148 + t168 * t150, t152 * t148 - t154 * t159, t160, t143 * t159 + t148 * t158, 0; 0, t151, 0, qJD(2) * t143, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:17:08
	% EndTime: 2019-12-29 18:17:08
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (260->40), mult. (402->58), div. (0->0), fcn. (344->8), ass. (0->34)
	t70 = qJ(2) + pkin(8);
	t68 = sin(t70);
	t105 = pkin(3) + pkin(4);
	t69 = cos(t70);
	t81 = t105 * t68 + sin(qJ(2)) * pkin(2) - qJ(4) * t69;
	t78 = -t81 * qJD(2) + t68 * qJD(4);
	t111 = t78 - (pkin(7) + r_i_i_C(3) - qJ(3) - pkin(6)) * qJD(1);
	t72 = sin(qJ(5));
	t75 = cos(qJ(5));
	t110 = -t68 * t75 + t69 * t72;
	t107 = qJD(2) - qJD(5);
	t86 = t68 * t72 + t69 * t75;
	t62 = t107 * t86;
	t109 = -qJ(4) * t68 - t105 * t69 - cos(qJ(2)) * pkin(2);
	t99 = qJD(2) * t68;
	t106 = t110 * qJD(5) + t75 * t99;
	t74 = sin(qJ(1));
	t101 = qJD(1) * t74;
	t77 = cos(qJ(1));
	t100 = qJD(1) * t77;
	t98 = qJD(2) * t69;
	t92 = t77 * t98;
	t84 = qJD(1) * t110;
	t57 = t62 * t77 + t74 * t84;
	t83 = qJD(1) * t86;
	t58 = t106 * t77 - t72 * t92 + t74 * t83;
	t91 = t57 * r_i_i_C(1) + t58 * r_i_i_C(2);
	t59 = -t74 * t62 + t77 * t84;
	t60 = t77 * t83 + (t72 * t98 - t106) * t74;
	t90 = -t59 * r_i_i_C(1) - t60 * r_i_i_C(2);
	t89 = -r_i_i_C(1) * t107 * t110 - t62 * r_i_i_C(2);
	t80 = qJD(3) + (-pkin(1) + t109) * qJD(1);
	t79 = t109 * qJD(2) + qJD(4) * t69;
	t1 = [-t60 * r_i_i_C(1) + t59 * r_i_i_C(2) - t111 * t74 + t80 * t77, t81 * t101 + t79 * t77 - t91, t100, -t68 * t101 + t92, t91; -t58 * r_i_i_C(1) + t57 * r_i_i_C(2) + t111 * t77 + t80 * t74, -t100 * t81 + t79 * t74 - t90, t101, t68 * t100 + t74 * t98, t90; 0, t78 - t89, 0, t99, t89;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end