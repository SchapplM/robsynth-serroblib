% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR14_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (41->8), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t69 = cos(pkin(5));
	t70 = sin(qJ(3));
	t73 = t69 * t70;
	t71 = cos(qJ(3));
	t72 = t69 * t71;
	t68 = sin(pkin(5));
	t67 = qJ(1) + qJ(2);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = t66 * t68;
	t63 = t65 * t68;
	t62 = -t65 * t73 + t66 * t71;
	t61 = -t65 * t72 - t66 * t70;
	t60 = -t65 * t71 - t66 * t73;
	t59 = t65 * t70 - t66 * t72;
	t1 = [t60, t60, t61, 0, 0; t62, t62, -t59, 0, 0; 0, 0, t68 * t71, 0, 0; t59, t59, -t62, 0, 0; t61, t61, t60, 0, 0; 0, 0, -t68 * t70, 0, 0; t64, t64, 0, 0, 0; t63, t63, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (168->16), mult. (76->16), div. (0->0), fcn. (96->9), ass. (0->24)
	t111 = qJ(3) + qJ(4);
	t106 = pkin(5) + t111;
	t104 = cos(t106) / 0.2e1;
	t116 = pkin(5) - t111;
	t105 = cos(t116);
	t118 = t105 / 0.2e1 + t104;
	t117 = sin(t106) / 0.2e1;
	t115 = sin(t116);
	t112 = qJ(1) + qJ(2);
	t108 = sin(t112);
	t109 = cos(t111);
	t110 = cos(t112);
	t98 = t117 - t115 / 0.2e1;
	t95 = t108 * t98 - t110 * t109;
	t114 = -t108 * t109 - t110 * t98;
	t113 = sin(pkin(5));
	t107 = sin(t111);
	t101 = t110 * t113;
	t100 = t108 * t113;
	t99 = t104 - t105 / 0.2e1;
	t97 = t117 + t115 / 0.2e1;
	t93 = -t110 * t107 - t108 * t118;
	t90 = t108 * t107 - t110 * t118;
	t1 = [t114, t114, t93, t93, 0; -t95, -t95, -t90, -t90, 0; 0, 0, t97, t97, 0; t90, t90, t95, t95, 0; t93, t93, t114, t114, 0; 0, 0, t99, t99, 0; t101, t101, 0, 0, 0; t100, t100, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (286->17), mult. (96->16), div. (0->0), fcn. (120->9), ass. (0->24)
	t119 = qJ(3) + qJ(4) + qJ(5);
	t114 = pkin(5) + t119;
	t128 = cos(t114) / 0.2e1;
	t127 = sin(t114) / 0.2e1;
	t126 = pkin(5) - t119;
	t123 = sin(t126);
	t108 = t127 - t123 / 0.2e1;
	t116 = cos(t119);
	t120 = qJ(1) + qJ(2);
	t117 = sin(t120);
	t118 = cos(t120);
	t125 = -t108 * t118 - t116 * t117;
	t105 = t108 * t117 - t116 * t118;
	t124 = cos(t126);
	t122 = t124 / 0.2e1 + t128;
	t121 = sin(pkin(5));
	t115 = sin(t119);
	t113 = t118 * t121;
	t112 = t117 * t121;
	t109 = t128 - t124 / 0.2e1;
	t107 = t127 + t123 / 0.2e1;
	t103 = -t118 * t115 - t117 * t122;
	t100 = t115 * t117 - t118 * t122;
	t1 = [t125, t125, t103, t103, t103; -t105, -t105, -t100, -t100, -t100; 0, 0, t107, t107, t107; t100, t100, t105, t105, t105; t103, t103, t125, t125, t125; 0, 0, t109, t109, t109; t113, t113, 0, 0, 0; t112, t112, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end