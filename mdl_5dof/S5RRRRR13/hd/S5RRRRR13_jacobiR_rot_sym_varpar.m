% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR13
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
%   Siehe auch: S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR13_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
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
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
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
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(1) + qJ(2) + qJ(3);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [-t19, -t19, -t19, 0, 0; t20, t20, t20, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, -t20, 0, 0; -t19, -t19, -t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (93->8), mult. (56->14), div. (0->0), fcn. (96->6), ass. (0->16)
	t74 = cos(pkin(5));
	t75 = sin(qJ(4));
	t78 = t74 * t75;
	t76 = cos(qJ(4));
	t77 = t74 * t76;
	t73 = sin(pkin(5));
	t72 = qJ(1) + qJ(2) + qJ(3);
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = t71 * t73;
	t68 = t70 * t73;
	t67 = -t70 * t78 + t71 * t76;
	t66 = -t70 * t77 - t71 * t75;
	t65 = -t70 * t76 - t71 * t78;
	t64 = t70 * t75 - t71 * t77;
	t1 = [t65, t65, t65, t66, 0; t67, t67, t67, -t64, 0; 0, 0, 0, t73 * t76, 0; t64, t64, t64, -t67, 0; t66, t66, t66, t65, 0; 0, 0, 0, -t73 * t75, 0; t69, t69, t69, 0, 0; t68, t68, t68, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (252->17), mult. (94->16), div. (0->0), fcn. (120->9), ass. (0->24)
	t117 = qJ(4) + qJ(5);
	t113 = pkin(5) + t117;
	t109 = cos(t113) / 0.2e1;
	t121 = pkin(5) - t117;
	t110 = cos(t121);
	t123 = t110 / 0.2e1 + t109;
	t122 = sin(t113) / 0.2e1;
	t120 = sin(t121);
	t103 = t122 - t120 / 0.2e1;
	t116 = qJ(1) + qJ(2) + qJ(3);
	t111 = sin(t116);
	t112 = cos(t116);
	t115 = cos(t117);
	t119 = -t112 * t103 - t111 * t115;
	t100 = t111 * t103 - t112 * t115;
	t118 = sin(pkin(5));
	t114 = sin(t117);
	t106 = t112 * t118;
	t105 = t111 * t118;
	t104 = t109 - t110 / 0.2e1;
	t102 = t122 + t120 / 0.2e1;
	t98 = -t111 * t123 - t112 * t114;
	t95 = t111 * t114 - t112 * t123;
	t1 = [t119, t119, t119, t98, t98; -t100, -t100, -t100, -t95, -t95; 0, 0, 0, t102, t102; t95, t95, t95, t100, t100; t98, t98, t98, t119, t119; 0, 0, 0, t104, t104; t106, t106, t106, 0, 0; t105, t105, t105, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end