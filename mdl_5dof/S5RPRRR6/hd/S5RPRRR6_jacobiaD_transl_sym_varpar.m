% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:46
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t68 = sin(qJ(3));
	t69 = cos(qJ(3));
	t70 = qJD(3) * (r_i_i_C(1) * t68 + r_i_i_C(2) * t69);
	t77 = -pkin(6) - r_i_i_C(3);
	t76 = qJD(1) * t68;
	t75 = qJD(1) * t69;
	t74 = qJD(3) * t68;
	t73 = qJD(3) * t69;
	t71 = r_i_i_C(1) * t69 - r_i_i_C(2) * t68 + pkin(2);
	t67 = qJ(1) + pkin(9);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [0, 0, -t70, 0, 0; t66 * t70 + (sin(qJ(1)) * pkin(1) + t77 * t66 + t71 * t65) * qJD(1), 0, (-t65 * t74 + t66 * t75) * r_i_i_C(2) + (t65 * t73 + t66 * t76) * r_i_i_C(1), 0, 0; t65 * t70 + (-cos(qJ(1)) * pkin(1) + t77 * t65 - t71 * t66) * qJD(1), 0, (t65 * t75 + t66 * t74) * r_i_i_C(2) + (t65 * t76 - t66 * t73) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (119->23), mult. (118->32), div. (0->0), fcn. (75->8), ass. (0->22)
	t104 = qJ(3) + qJ(4);
	t100 = sin(t104);
	t101 = cos(t104);
	t129 = r_i_i_C(1) * t100 + r_i_i_C(2) * t101;
	t128 = -r_i_i_C(1) * t101 + r_i_i_C(2) * t100;
	t105 = sin(qJ(3));
	t102 = qJD(3) + qJD(4);
	t111 = t129 * t102;
	t108 = qJD(3) * t105 * pkin(3) + t111;
	t127 = t129 * qJD(1);
	t126 = t128 * t102;
	t121 = -r_i_i_C(3) - pkin(7) - pkin(6);
	t120 = qJD(1) * t105;
	t106 = cos(qJ(3));
	t119 = qJD(3) * t106;
	t112 = pkin(3) * t106 + pkin(2) - t128;
	t103 = qJ(1) + pkin(9);
	t98 = sin(t103);
	t99 = cos(t103);
	t110 = t126 * t99 + t127 * t98;
	t109 = -t126 * t98 + t127 * t99;
	t1 = [0, 0, -t108, -t111, 0; t108 * t99 + (sin(qJ(1)) * pkin(1) + t121 * t99 + t112 * t98) * qJD(1), 0, (t119 * t98 + t120 * t99) * pkin(3) + t109, t109, 0; t108 * t98 + (-cos(qJ(1)) * pkin(1) + t121 * t98 - t112 * t99) * qJD(1), 0, (-t119 * t99 + t120 * t98) * pkin(3) + t110, t110, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:15
	% EndTime: 2019-10-24 10:46:15
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (252->36), mult. (174->46), div. (0->0), fcn. (111->10), ass. (0->35)
	t116 = qJ(3) + qJ(4);
	t112 = qJ(5) + t116;
	t105 = sin(t112);
	t106 = cos(t112);
	t146 = r_i_i_C(1) * t105 + r_i_i_C(2) * t106;
	t117 = sin(qJ(3));
	t113 = qJD(3) + qJD(4);
	t109 = qJD(5) + t113;
	t123 = t146 * t109;
	t110 = sin(t116);
	t143 = pkin(4) * t110;
	t119 = -t113 * t143 - t123;
	t137 = pkin(3) * qJD(3);
	t145 = -t117 * t137 + t119;
	t144 = t146 * qJD(1);
	t141 = r_i_i_C(1) * t106;
	t140 = r_i_i_C(2) * t105;
	t138 = -r_i_i_C(3) - pkin(8) - pkin(7) - pkin(6);
	t111 = cos(t116);
	t136 = t111 * t113;
	t135 = qJD(1) * (-t117 * pkin(3) - t143);
	t134 = qJD(1) * t110;
	t132 = t109 * t141;
	t131 = t109 * t140;
	t114 = qJ(1) + pkin(9);
	t107 = sin(t114);
	t108 = cos(t114);
	t130 = t144 * t107 + t108 * t131;
	t126 = t107 * t132 + t144 * t108;
	t118 = cos(qJ(3));
	t124 = t118 * pkin(3) + pkin(4) * t111 + pkin(2) - t140 + t141;
	t122 = -t108 * t132 + t130;
	t121 = -t107 * t131 + t126;
	t98 = -pkin(4) * t136 - t118 * t137;
	t1 = [0, 0, t145, t119, -t123; -t145 * t108 + (sin(qJ(1)) * pkin(1) + t138 * t108 + t124 * t107) * qJD(1), 0, -t108 * t135 + (-t98 - t131) * t107 + t126, (t107 * t136 + t108 * t134) * pkin(4) + t121, t121; -t145 * t107 + (-cos(qJ(1)) * pkin(1) + t138 * t107 - t124 * t108) * qJD(1), 0, -t107 * t135 + (t98 - t132) * t108 + t130, (t107 * t134 - t108 * t136) * pkin(4) + t122, t122;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end