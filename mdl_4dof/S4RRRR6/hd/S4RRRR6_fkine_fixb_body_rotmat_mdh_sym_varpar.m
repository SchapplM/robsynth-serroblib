% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:41
	% EndTime: 2020-11-04 19:51:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:41
	% EndTime: 2020-11-04 19:51:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:41
	% EndTime: 2020-11-04 19:51:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t65 = sin(pkin(4));
	t68 = sin(qJ(1));
	t76 = t68 * t65;
	t67 = sin(qJ(2));
	t75 = t68 * t67;
	t69 = cos(qJ(2));
	t74 = t68 * t69;
	t70 = cos(qJ(1));
	t73 = t70 * t65;
	t72 = t70 * t67;
	t71 = t70 * t69;
	t66 = cos(pkin(4));
	t1 = [-t66 * t75 + t71, -t66 * t74 - t72, t76, t70 * pkin(1) + pkin(6) * t76 + 0; t66 * t72 + t74, t66 * t71 - t75, -t73, t68 * pkin(1) - pkin(6) * t73 + 0; t65 * t67, t65 * t69, t66, t66 * pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:41
	% EndTime: 2020-11-04 19:51:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t80 = sin(pkin(4));
	t85 = cos(qJ(3));
	t95 = t80 * t85;
	t82 = sin(qJ(3));
	t94 = t82 * t80;
	t83 = sin(qJ(2));
	t93 = t83 * t85;
	t84 = sin(qJ(1));
	t92 = t84 * t83;
	t86 = cos(qJ(2));
	t91 = t84 * t86;
	t87 = cos(qJ(1));
	t90 = t87 * t83;
	t89 = t87 * t86;
	t88 = pkin(2) * t83 - pkin(7) * t86;
	t81 = cos(pkin(4));
	t79 = t86 * pkin(2) + t83 * pkin(7) + pkin(1);
	t78 = t81 * t90 + t91;
	t77 = t80 * pkin(6) - t88 * t81;
	t1 = [(-t81 * t93 + t94) * t84 + t85 * t89, (t81 * t92 - t89) * t82 + t84 * t95, t81 * t91 + t90, t77 * t84 + t79 * t87 + 0; t78 * t85 - t87 * t94, -t78 * t82 - t87 * t95, -t81 * t89 + t92, -t77 * t87 + t79 * t84 + 0; t80 * t93 + t81 * t82, t81 * t85 - t83 * t94, -t80 * t86, t81 * pkin(6) + t88 * t80 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:41
	% EndTime: 2020-11-04 19:51:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t103 = sin(pkin(4));
	t110 = cos(qJ(3));
	t124 = t103 * t110;
	t105 = sin(qJ(4));
	t111 = cos(qJ(2));
	t123 = t105 * t111;
	t106 = sin(qJ(3));
	t122 = t106 * t103;
	t107 = sin(qJ(2));
	t121 = t107 * t110;
	t108 = sin(qJ(1));
	t120 = t108 * t111;
	t109 = cos(qJ(4));
	t119 = t109 * t111;
	t118 = t110 * t111;
	t112 = cos(qJ(1));
	t117 = t112 * t107;
	t116 = t112 * t111;
	t101 = t110 * pkin(3) + t106 * pkin(8) + pkin(2);
	t115 = t111 * pkin(7) - t101 * t107;
	t114 = t106 * pkin(3) - t110 * pkin(8) + pkin(6);
	t104 = cos(pkin(4));
	t99 = t104 * t121 - t122;
	t113 = t104 * t119 + t99 * t105;
	t100 = t105 * t118 - t109 * t107;
	t98 = t103 * t121 + t104 * t106;
	t97 = t107 * pkin(7) + t101 * t111 + pkin(1);
	t96 = t103 * t114 + t115 * t104;
	t1 = [(-t99 * t108 + t110 * t116) * t109 + (t104 * t120 + t117) * t105, -t112 * t100 + t113 * t108, (-t108 * t104 * t107 + t116) * t106 - t108 * t124, t96 * t108 + t97 * t112 + 0; (-t104 * t123 + t99 * t109) * t112 + t108 * (t107 * t105 + t109 * t118), -t108 * t100 - t113 * t112, (t104 * t117 + t120) * t106 + t112 * t124, t97 * t108 - t96 * t112 + 0; -t103 * t123 + t98 * t109, -t103 * t119 - t98 * t105, -t104 * t110 + t107 * t122, -t115 * t103 + t114 * t104 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end