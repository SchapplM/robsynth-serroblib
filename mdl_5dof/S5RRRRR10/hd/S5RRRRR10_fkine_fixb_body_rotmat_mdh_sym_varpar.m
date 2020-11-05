% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t1 = [t72, -t71, 0, 0; t71, t72, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t73 = sin(pkin(5));
	t76 = sin(qJ(1));
	t84 = t76 * t73;
	t75 = sin(qJ(2));
	t83 = t76 * t75;
	t77 = cos(qJ(2));
	t82 = t76 * t77;
	t78 = cos(qJ(1));
	t81 = t78 * t73;
	t80 = t78 * t75;
	t79 = t78 * t77;
	t74 = cos(pkin(5));
	t1 = [-t74 * t83 + t79, -t74 * t82 - t80, t84, t78 * pkin(1) + pkin(7) * t84 + 0; t74 * t80 + t82, t74 * t79 - t83, -t81, t76 * pkin(1) - pkin(7) * t81 + 0; t73 * t75, t73 * t77, t74, t74 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t88 = sin(pkin(5));
	t93 = cos(qJ(3));
	t103 = t88 * t93;
	t90 = sin(qJ(3));
	t102 = t90 * t88;
	t91 = sin(qJ(2));
	t101 = t91 * t93;
	t92 = sin(qJ(1));
	t100 = t92 * t91;
	t94 = cos(qJ(2));
	t99 = t92 * t94;
	t95 = cos(qJ(1));
	t98 = t95 * t91;
	t97 = t95 * t94;
	t96 = pkin(2) * t91 - pkin(8) * t94;
	t89 = cos(pkin(5));
	t87 = t94 * pkin(2) + t91 * pkin(8) + pkin(1);
	t86 = -t89 * t98 - t99;
	t85 = t88 * pkin(7) - t96 * t89;
	t1 = [(-t89 * t101 + t102) * t92 + t93 * t97, (t89 * t100 - t97) * t90 + t92 * t103, t89 * t99 + t98, t85 * t92 + t87 * t95 + 0; -t95 * t102 - t86 * t93, -t95 * t103 + t86 * t90, -t89 * t97 + t100, -t85 * t95 + t87 * t92 + 0; t88 * t101 + t89 * t90, -t91 * t102 + t89 * t93, -t88 * t94, t89 * pkin(7) + t96 * t88 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t112 = sin(pkin(5));
	t115 = sin(qJ(2));
	t128 = t112 * t115;
	t116 = sin(qJ(1));
	t127 = t112 * t116;
	t118 = cos(qJ(1));
	t126 = t112 * t118;
	t125 = t116 * t115;
	t117 = cos(qJ(2));
	t124 = t116 * t117;
	t123 = t118 * t115;
	t122 = t118 * t117;
	t121 = sin(qJ(3)) * pkin(3) + pkin(7);
	t108 = cos(qJ(3)) * pkin(3) + pkin(2);
	t119 = pkin(9) + pkin(8);
	t120 = t108 * t115 - t119 * t117;
	t113 = cos(pkin(5));
	t111 = qJ(3) + qJ(4);
	t110 = cos(t111);
	t109 = sin(t111);
	t107 = t113 * t123 + t124;
	t106 = t113 * t125 - t122;
	t105 = t108 * t117 + t119 * t115 + pkin(1);
	t104 = t112 * t121 - t120 * t113;
	t1 = [-t106 * t110 + t109 * t127, t106 * t109 + t110 * t127, t113 * t124 + t123, t104 * t116 + t105 * t118 + 0; t107 * t110 - t109 * t126, -t107 * t109 - t110 * t126, -t113 * t122 + t125, -t104 * t118 + t105 * t116 + 0; t113 * t109 + t110 * t128, -t109 * t128 + t113 * t110, -t112 * t117, t120 * t112 + t121 * t113 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:15
	% EndTime: 2020-11-04 20:51:15
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (100->38), mult. (139->56), div. (0->0), fcn. (183->12), ass. (0->36)
	t145 = sin(pkin(5));
	t149 = sin(qJ(2));
	t164 = t145 * t149;
	t150 = sin(qJ(1));
	t163 = t145 * t150;
	t152 = cos(qJ(2));
	t162 = t145 * t152;
	t153 = cos(qJ(1));
	t161 = t145 * t153;
	t160 = t150 * t149;
	t159 = t150 * t152;
	t158 = t153 * t149;
	t157 = t153 * t152;
	t156 = sin(qJ(3)) * pkin(3) + pkin(7);
	t141 = cos(qJ(3)) * pkin(3) + pkin(2);
	t154 = pkin(9) + pkin(8);
	t155 = t141 * t149 - t154 * t152;
	t151 = cos(qJ(5));
	t147 = sin(qJ(5));
	t146 = cos(pkin(5));
	t144 = qJ(3) + qJ(4);
	t143 = cos(t144);
	t142 = sin(t144);
	t140 = t146 * t159 + t158;
	t139 = t146 * t158 + t159;
	t138 = -t146 * t157 + t160;
	t137 = t146 * t160 - t157;
	t136 = t141 * t152 + t154 * t149 + pkin(1);
	t135 = t146 * t142 + t143 * t164;
	t134 = t142 * t164 - t146 * t143;
	t133 = t145 * t156 - t146 * t155;
	t132 = -t137 * t143 + t142 * t163;
	t131 = t139 * t143 - t142 * t161;
	t130 = t139 * t142 + t143 * t161;
	t129 = t137 * t142 + t143 * t163;
	t1 = [t132 * t151 + t140 * t147, -t132 * t147 + t140 * t151, -t129, t132 * pkin(4) - t129 * pkin(10) + t133 * t150 + t136 * t153 + 0; t131 * t151 + t138 * t147, -t131 * t147 + t138 * t151, t130, t131 * pkin(4) + t130 * pkin(10) - t133 * t153 + t136 * t150 + 0; t135 * t151 - t147 * t162, -t135 * t147 - t151 * t162, t134, t135 * pkin(4) + t134 * pkin(10) + t155 * t145 + t146 * t156 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end