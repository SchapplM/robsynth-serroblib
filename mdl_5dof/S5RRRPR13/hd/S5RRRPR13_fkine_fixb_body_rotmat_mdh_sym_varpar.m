% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t77 = cos(qJ(1));
	t76 = sin(qJ(1));
	t1 = [t77, -t76, 0, 0; t76, t77, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t78 = sin(pkin(5));
	t81 = sin(qJ(1));
	t89 = t81 * t78;
	t80 = sin(qJ(2));
	t88 = t81 * t80;
	t82 = cos(qJ(2));
	t87 = t81 * t82;
	t83 = cos(qJ(1));
	t86 = t83 * t78;
	t85 = t83 * t80;
	t84 = t83 * t82;
	t79 = cos(pkin(5));
	t1 = [-t79 * t88 + t84, -t79 * t87 - t85, t89, t83 * pkin(1) + pkin(7) * t89 + 0; t79 * t85 + t87, t79 * t84 - t88, -t86, t81 * pkin(1) - pkin(7) * t86 + 0; t78 * t80, t78 * t82, t79, t79 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->38), div. (0->0), fcn. (81->8), ass. (0->20)
	t93 = sin(pkin(5));
	t95 = sin(qJ(3));
	t108 = t95 * t93;
	t96 = sin(qJ(2));
	t98 = cos(qJ(3));
	t107 = t96 * t98;
	t97 = sin(qJ(1));
	t106 = t97 * t96;
	t99 = cos(qJ(2));
	t105 = t97 * t99;
	t100 = cos(qJ(1));
	t104 = t100 * t93;
	t103 = t100 * t96;
	t102 = t100 * t99;
	t101 = pkin(2) * t96 - pkin(8) * t99;
	t94 = cos(pkin(5));
	t92 = t99 * pkin(2) + t96 * pkin(8) + pkin(1);
	t91 = -t94 * t103 - t105;
	t90 = t93 * pkin(7) - t101 * t94;
	t1 = [(-t94 * t107 + t108) * t97 + t98 * t102, (t94 * t106 - t102) * t95 + t97 * t93 * t98, t94 * t105 + t103, t92 * t100 + t90 * t97 + 0; -t95 * t104 - t91 * t98, -t98 * t104 + t91 * t95, -t94 * t102 + t106, -t90 * t100 + t92 * t97 + 0; t93 * t107 + t94 * t95, -t96 * t108 + t94 * t98, -t93 * t99, t94 * pkin(7) + t101 * t93 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t113 = sin(pkin(5));
	t118 = cos(qJ(3));
	t129 = t113 * t118;
	t115 = sin(qJ(3));
	t128 = t115 * t113;
	t116 = sin(qJ(2));
	t127 = t116 * t118;
	t117 = sin(qJ(1));
	t126 = t117 * t116;
	t119 = cos(qJ(2));
	t125 = t117 * t119;
	t120 = cos(qJ(1));
	t124 = t120 * t116;
	t123 = t120 * t119;
	t112 = t118 * pkin(3) + qJ(4) * t115 + pkin(2);
	t122 = t119 * pkin(8) - t112 * t116;
	t121 = t115 * pkin(3) - qJ(4) * t118 + pkin(7);
	t114 = cos(pkin(5));
	t111 = -t114 * t124 - t125;
	t110 = t116 * pkin(8) + t112 * t119 + pkin(1);
	t109 = t113 * t121 + t122 * t114;
	t1 = [t114 * t125 + t124, (t114 * t127 - t128) * t117 - t118 * t123, (-t114 * t126 + t123) * t115 - t117 * t129, t109 * t117 + t110 * t120 + 0; -t114 * t123 + t126, t111 * t118 + t120 * t128, -t111 * t115 + t120 * t129, -t109 * t120 + t110 * t117 + 0; -t113 * t119, -t113 * t127 - t114 * t115, -t114 * t118 + t116 * t128, -t122 * t113 + t121 * t114 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:46
	% EndTime: 2020-11-04 20:44:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (70->36), mult. (112->58), div. (0->0), fcn. (146->10), ass. (0->31)
	t137 = cos(pkin(5));
	t143 = cos(qJ(3));
	t159 = t137 * t143;
	t138 = sin(qJ(5));
	t144 = cos(qJ(2));
	t158 = t138 * t144;
	t136 = sin(pkin(5));
	t139 = sin(qJ(3));
	t157 = t139 * t136;
	t140 = sin(qJ(2));
	t156 = t139 * t140;
	t141 = sin(qJ(1));
	t155 = t141 * t144;
	t142 = cos(qJ(5));
	t154 = t142 * t144;
	t153 = t143 * t136;
	t145 = cos(qJ(1));
	t152 = t144 * t145;
	t151 = t145 * t140;
	t147 = pkin(3) + pkin(9);
	t135 = qJ(4) * t139 + t147 * t143 + pkin(2);
	t146 = pkin(4) + pkin(8);
	t150 = t135 * t140 - t146 * t144;
	t149 = qJ(4) * t143 - t147 * t139 - pkin(7);
	t132 = t137 * t156 + t153;
	t148 = -t132 * t138 + t137 * t154;
	t134 = t139 * t158 + t142 * t140;
	t133 = t136 * t156 - t159;
	t131 = t135 * t144 + t146 * t140 + pkin(1);
	t130 = -t136 * t149 - t150 * t137;
	t1 = [t145 * t134 + t148 * t141, (-t132 * t141 + t139 * t152) * t142 - (t137 * t155 + t151) * t138, (-t140 * t159 + t157) * t141 + t143 * t152, t130 * t141 + t131 * t145 + 0; t141 * t134 - t148 * t145, (t132 * t142 + t137 * t158) * t145 + t141 * (-t138 * t140 + t139 * t154), (t137 * t151 + t155) * t143 - t145 * t157, -t130 * t145 + t131 * t141 + 0; t133 * t138 - t136 * t154, t133 * t142 + t136 * t158, t137 * t139 + t140 * t153, t150 * t136 - t149 * t137 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end