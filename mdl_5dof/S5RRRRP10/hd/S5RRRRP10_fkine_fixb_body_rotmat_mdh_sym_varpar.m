% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t79 = cos(qJ(1));
	t78 = sin(qJ(1));
	t1 = [t79, -t78, 0, 0; t78, t79, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t80 = sin(pkin(5));
	t83 = sin(qJ(1));
	t91 = t83 * t80;
	t82 = sin(qJ(2));
	t90 = t83 * t82;
	t84 = cos(qJ(2));
	t89 = t83 * t84;
	t85 = cos(qJ(1));
	t88 = t85 * t80;
	t87 = t85 * t82;
	t86 = t85 * t84;
	t81 = cos(pkin(5));
	t1 = [-t81 * t90 + t86, -t81 * t89 - t87, t91, t85 * pkin(1) + pkin(7) * t91 + 0; t81 * t87 + t89, t81 * t86 - t90, -t88, t83 * pkin(1) - pkin(7) * t88 + 0; t80 * t82, t80 * t84, t81, t81 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->20)
	t95 = sin(pkin(5));
	t97 = sin(qJ(3));
	t110 = t97 * t95;
	t98 = sin(qJ(2));
	t99 = sin(qJ(1));
	t109 = t99 * t98;
	t102 = cos(qJ(1));
	t108 = t102 * t98;
	t100 = cos(qJ(3));
	t107 = t95 * t100;
	t96 = cos(pkin(5));
	t106 = t96 * t100;
	t101 = cos(qJ(2));
	t105 = t99 * t101;
	t104 = t102 * t101;
	t103 = pkin(2) * t98 - pkin(8) * t101;
	t94 = t101 * pkin(2) + t98 * pkin(8) + pkin(1);
	t93 = -t96 * t108 - t105;
	t92 = t95 * pkin(7) - t103 * t96;
	t1 = [(-t98 * t106 + t110) * t99 + t100 * t104, (t96 * t109 - t104) * t97 + t99 * t107, t96 * t105 + t108, t94 * t102 + t92 * t99 + 0; -t93 * t100 - t102 * t110, -t102 * t107 + t93 * t97, -t96 * t104 + t109, -t92 * t102 + t94 * t99 + 0; t98 * t107 + t96 * t97, -t98 * t110 + t106, -t95 * t101, t96 * pkin(7) + t103 * t95 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t118 = sin(pkin(5));
	t125 = cos(qJ(3));
	t139 = t118 * t125;
	t120 = sin(qJ(4));
	t126 = cos(qJ(2));
	t138 = t120 * t126;
	t121 = sin(qJ(3));
	t137 = t121 * t118;
	t122 = sin(qJ(2));
	t136 = t122 * t125;
	t123 = sin(qJ(1));
	t135 = t123 * t126;
	t124 = cos(qJ(4));
	t134 = t124 * t126;
	t133 = t125 * t126;
	t127 = cos(qJ(1));
	t132 = t127 * t122;
	t131 = t127 * t126;
	t116 = t125 * pkin(3) + t121 * pkin(9) + pkin(2);
	t130 = t126 * pkin(8) - t116 * t122;
	t129 = t121 * pkin(3) - t125 * pkin(9) + pkin(7);
	t119 = cos(pkin(5));
	t114 = t119 * t136 - t137;
	t128 = t114 * t120 + t119 * t134;
	t115 = t120 * t133 - t124 * t122;
	t113 = t118 * t136 + t119 * t121;
	t112 = t122 * pkin(8) + t116 * t126 + pkin(1);
	t111 = t118 * t129 + t130 * t119;
	t1 = [(-t114 * t123 + t125 * t131) * t124 + (t119 * t135 + t132) * t120, -t127 * t115 + t128 * t123, (-t123 * t119 * t122 + t131) * t121 - t123 * t139, t111 * t123 + t112 * t127 + 0; (t114 * t124 - t119 * t138) * t127 + t123 * (t120 * t122 + t124 * t133), -t123 * t115 - t128 * t127, (t119 * t132 + t135) * t121 + t127 * t139, -t111 * t127 + t112 * t123 + 0; t113 * t124 - t118 * t138, -t113 * t120 - t118 * t134, -t119 * t125 + t122 * t137, -t130 * t118 + t129 * t119 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:15
	% EndTime: 2020-11-04 20:47:15
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (78->37), mult. (125->62), div. (0->0), fcn. (159->10), ass. (0->32)
	t149 = sin(pkin(5));
	t157 = cos(qJ(3));
	t171 = t149 * t157;
	t152 = sin(qJ(4));
	t158 = cos(qJ(2));
	t170 = t152 * t158;
	t153 = sin(qJ(3));
	t169 = t153 * t149;
	t154 = sin(qJ(2));
	t168 = t154 * t157;
	t155 = sin(qJ(1));
	t167 = t155 * t158;
	t156 = cos(qJ(4));
	t166 = t156 * t158;
	t165 = t157 * t158;
	t159 = cos(qJ(1));
	t164 = t159 * t154;
	t163 = t159 * t158;
	t148 = t156 * pkin(4) + pkin(3);
	t151 = qJ(5) + pkin(9);
	t142 = t148 * t157 + t151 * t153 + pkin(2);
	t147 = t152 * pkin(4) + pkin(8);
	t162 = t142 * t154 - t147 * t158;
	t161 = t148 * t153 - t151 * t157 + pkin(7);
	t150 = cos(pkin(5));
	t144 = t150 * t168 - t169;
	t160 = t144 * t152 + t150 * t166;
	t145 = t152 * t165 - t156 * t154;
	t143 = t149 * t168 + t150 * t153;
	t141 = t142 * t158 + t147 * t154 + pkin(1);
	t140 = t149 * t161 - t162 * t150;
	t1 = [(-t144 * t155 + t157 * t163) * t156 + (t150 * t167 + t164) * t152, -t159 * t145 + t160 * t155, (-t155 * t150 * t154 + t163) * t153 - t155 * t171, t140 * t155 + t141 * t159 + 0; (t144 * t156 - t150 * t170) * t159 + t155 * (t152 * t154 + t156 * t165), -t155 * t145 - t160 * t159, (t150 * t164 + t167) * t153 + t159 * t171, -t140 * t159 + t141 * t155 + 0; t143 * t156 - t149 * t170, -t143 * t152 - t149 * t166, -t150 * t157 + t154 * t169, t162 * t149 + t161 * t150 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end