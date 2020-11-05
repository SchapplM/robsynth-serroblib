% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t1 = [t80, -t79, 0, 0; t79, t80, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t81 = sin(pkin(5));
	t84 = sin(qJ(1));
	t92 = t84 * t81;
	t83 = sin(qJ(2));
	t91 = t84 * t83;
	t85 = cos(qJ(2));
	t90 = t84 * t85;
	t86 = cos(qJ(1));
	t89 = t86 * t81;
	t88 = t86 * t83;
	t87 = t86 * t85;
	t82 = cos(pkin(5));
	t1 = [-t82 * t91 + t87, -t82 * t90 - t88, t92, t86 * pkin(1) + pkin(7) * t92 + 0; t82 * t88 + t90, t82 * t87 - t91, -t89, t84 * pkin(1) - pkin(7) * t89 + 0; t81 * t83, t81 * t85, t82, t82 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->20)
	t96 = sin(pkin(5));
	t98 = sin(qJ(3));
	t111 = t98 * t96;
	t100 = sin(qJ(1));
	t99 = sin(qJ(2));
	t110 = t100 * t99;
	t103 = cos(qJ(1));
	t109 = t103 * t99;
	t101 = cos(qJ(3));
	t108 = t96 * t101;
	t97 = cos(pkin(5));
	t107 = t97 * t101;
	t102 = cos(qJ(2));
	t106 = t100 * t102;
	t105 = t103 * t102;
	t104 = pkin(2) * t99 - pkin(8) * t102;
	t95 = t102 * pkin(2) + t99 * pkin(8) + pkin(1);
	t94 = t97 * t109 + t106;
	t93 = t96 * pkin(7) - t104 * t97;
	t1 = [(-t99 * t107 + t111) * t100 + t101 * t105, (t97 * t110 - t105) * t98 + t100 * t108, t97 * t106 + t109, t93 * t100 + t95 * t103 + 0; t94 * t101 - t103 * t111, -t103 * t108 - t94 * t98, -t97 * t105 + t110, t95 * t100 - t93 * t103 + 0; t99 * t108 + t97 * t98, -t99 * t111 + t107, -t96 * t102, t97 * pkin(7) + t104 * t96 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t120 = sin(pkin(5));
	t124 = sin(qJ(2));
	t136 = t120 * t124;
	t125 = sin(qJ(1));
	t135 = t120 * t125;
	t127 = cos(qJ(1));
	t134 = t120 * t127;
	t133 = t125 * t124;
	t126 = cos(qJ(2));
	t132 = t125 * t126;
	t131 = t127 * t124;
	t130 = t127 * t126;
	t129 = sin(qJ(3)) * pkin(3) + pkin(7);
	t116 = cos(qJ(3)) * pkin(3) + pkin(2);
	t122 = qJ(4) + pkin(8);
	t128 = t116 * t124 - t122 * t126;
	t121 = cos(pkin(5));
	t119 = qJ(3) + pkin(10);
	t118 = cos(t119);
	t117 = sin(t119);
	t115 = -t121 * t133 + t130;
	t114 = t121 * t131 + t132;
	t113 = t116 * t126 + t122 * t124 + pkin(1);
	t112 = t120 * t129 - t128 * t121;
	t1 = [t115 * t118 + t117 * t135, -t115 * t117 + t118 * t135, t121 * t132 + t131, t112 * t125 + t113 * t127 + 0; t114 * t118 - t117 * t134, -t114 * t117 - t118 * t134, -t121 * t130 + t133, -t112 * t127 + t113 * t125 + 0; t121 * t117 + t118 * t136, -t117 * t136 + t121 * t118, -t120 * t126, t128 * t120 + t129 * t121 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:56
	% EndTime: 2020-11-04 20:43:56
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (104->39), mult. (152->61), div. (0->0), fcn. (186->14), ass. (0->40)
	t153 = sin(pkin(10));
	t155 = cos(pkin(10));
	t148 = t155 * pkin(4) + t153 * pkin(9) + pkin(3);
	t159 = sin(qJ(3));
	t145 = t148 * t159 + pkin(7);
	t149 = -t153 * pkin(4) + t155 * pkin(9);
	t146 = t149 * t159 + pkin(2);
	t154 = sin(pkin(5));
	t156 = cos(pkin(5));
	t160 = sin(qJ(2));
	t163 = cos(qJ(3));
	t157 = qJ(4) + pkin(8);
	t164 = cos(qJ(2));
	t174 = t157 * t164;
	t179 = (t146 * t160 - t174) * t156 + (t156 * t148 * t160 + t154 * t149) * t163 - t154 * t145;
	t178 = t154 * t160;
	t161 = sin(qJ(1));
	t177 = t154 * t161;
	t176 = t154 * t164;
	t165 = cos(qJ(1));
	t175 = t154 * t165;
	t173 = t161 * t160;
	t172 = t161 * t164;
	t171 = t165 * t160;
	t170 = t165 * t164;
	t142 = t156 * t171 + t172;
	t152 = qJ(3) + pkin(10);
	t150 = sin(t152);
	t151 = cos(t152);
	t167 = -t142 * t151 + t150 * t175;
	t144 = -t156 * t173 + t170;
	t166 = t144 * t151 + t150 * t177;
	t162 = cos(qJ(5));
	t158 = sin(qJ(5));
	t143 = t156 * t172 + t171;
	t141 = t156 * t170 - t173;
	t140 = t156 * t150 + t151 * t178;
	t139 = t148 * t163 + t146;
	t137 = t139 * t164 + t157 * t160 + pkin(1);
	t1 = [t143 * t158 + t166 * t162, t143 * t162 - t166 * t158, t144 * t150 - t151 * t177, t137 * t165 - t179 * t161 + 0; -t158 * t141 - t167 * t162, -t162 * t141 + t167 * t158, t142 * t150 + t151 * t175, t137 * t161 + t179 * t165 + 0; t140 * t162 - t158 * t176, -t140 * t158 - t162 * t176, t150 * t178 - t156 * t151, pkin(6) + 0 + (t139 * t160 - t174) * t154 + (-t149 * t163 + t145) * t156; 0, 0, 0, 1;];
	Tc_mdh = t1;
end