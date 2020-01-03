% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->17), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t76 = pkin(6) + r_i_i_C(3);
	t68 = sin(qJ(3));
	t75 = qJD(1) * t68;
	t69 = cos(qJ(3));
	t74 = qJD(1) * t69;
	t73 = qJD(3) * t68;
	t72 = qJD(3) * t69;
	t71 = r_i_i_C(1) * t69 - r_i_i_C(2) * t68 + pkin(2);
	t70 = (-r_i_i_C(1) * t68 - r_i_i_C(2) * t69) * qJD(3);
	t67 = qJ(1) + pkin(8);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [0, 0, t70, 0, 0; t66 * t70 + (-sin(qJ(1)) * pkin(1) + t76 * t66 - t71 * t65) * qJD(1), 0, (t65 * t73 - t66 * t74) * r_i_i_C(2) + (-t65 * t72 - t66 * t75) * r_i_i_C(1), 0, 0; t65 * t70 + (cos(qJ(1)) * pkin(1) + t76 * t65 + t71 * t66) * qJD(1), 0, (-t65 * t74 - t66 * t73) * r_i_i_C(2) + (-t65 * t75 + t66 * t72) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (119->28), mult. (118->43), div. (0->0), fcn. (75->8), ass. (0->26)
	t100 = qJ(3) + qJ(4);
	t96 = sin(t100);
	t98 = qJD(3) + qJD(4);
	t117 = t96 * t98;
	t97 = cos(t100);
	t116 = t97 * t98;
	t115 = r_i_i_C(3) + pkin(7) + pkin(6);
	t114 = pkin(3) * qJD(3);
	t113 = qJD(1) * t96;
	t112 = qJD(1) * t97;
	t111 = r_i_i_C(1) * t116;
	t110 = r_i_i_C(2) * t117;
	t102 = cos(qJ(3));
	t109 = t102 * t114;
	t108 = -r_i_i_C(1) * t96 - r_i_i_C(2) * t97;
	t107 = t102 * pkin(3) + r_i_i_C(1) * t97 - r_i_i_C(2) * t96 + pkin(2);
	t106 = t108 * t98;
	t101 = sin(qJ(3));
	t105 = qJD(1) * (-pkin(3) * t101 + t108);
	t104 = -t101 * t114 + t106;
	t99 = qJ(1) + pkin(8);
	t95 = cos(t99);
	t94 = sin(t99);
	t92 = t95 * t111;
	t91 = t94 * t110;
	t1 = [0, 0, t104, t106, 0; t104 * t95 + (-sin(qJ(1)) * pkin(1) + t115 * t95 - t107 * t94) * qJD(1), 0, t91 + (-t109 - t111) * t94 + t95 * t105, -t95 * r_i_i_C(2) * t112 + t91 + (-t95 * t113 - t94 * t116) * r_i_i_C(1), 0; t104 * t94 + (cos(qJ(1)) * pkin(1) + t115 * t94 + t107 * t95) * qJD(1), 0, t92 + (t109 - t110) * t95 + t94 * t105, -t94 * r_i_i_C(1) * t113 + t92 + (-t94 * t112 - t95 * t117) * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (168->36), mult. (150->48), div. (0->0), fcn. (97->8), ass. (0->30)
	t103 = qJ(3) + qJ(4);
	t98 = sin(t103);
	t99 = cos(t103);
	t124 = pkin(4) * t99 - r_i_i_C(2) * t98;
	t123 = -pkin(4) - r_i_i_C(1);
	t120 = r_i_i_C(2) * t99;
	t119 = r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6);
	t118 = pkin(3) * qJD(3);
	t101 = qJD(3) + qJD(4);
	t117 = t101 * t98;
	t116 = t101 * t99;
	t115 = r_i_i_C(1) * t116;
	t114 = r_i_i_C(2) * t117;
	t104 = sin(qJ(3));
	t113 = t104 * t118;
	t112 = -r_i_i_C(1) * t98 - t120;
	t105 = cos(qJ(3));
	t111 = t105 * pkin(3) + r_i_i_C(1) * t99 + pkin(2) + t124;
	t110 = t123 * t98 - t120;
	t109 = qJD(1) * (-t104 * pkin(3) - pkin(4) * t98 + t112);
	t108 = t110 * t101;
	t107 = qJD(1) * t110;
	t106 = -pkin(4) * t117 + t112 * t101 - t113;
	t102 = qJ(1) + pkin(8);
	t97 = cos(t102);
	t96 = sin(t102);
	t93 = t97 * t115;
	t92 = t96 * t114;
	t91 = -pkin(4) * t116 - t105 * t118;
	t1 = [0, 0, t108 - t113, t108, 0; t96 * qJD(5) + t106 * t97 + (-sin(qJ(1)) * pkin(1) + t119 * t97 - t111 * t96) * qJD(1), 0, t92 + (t91 - t115) * t96 + t97 * t109, t123 * t96 * t116 + t97 * t107 + t92, qJD(1) * t96; -t97 * qJD(5) + t106 * t96 + (cos(qJ(1)) * pkin(1) + t119 * t96 + t111 * t97) * qJD(1), 0, t93 + (-t91 - t114) * t97 + t96 * t109, t124 * t97 * t101 + t96 * t107 + t93, -qJD(1) * t97;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end