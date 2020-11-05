% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t94 = cos(qJ(1));
	t93 = sin(qJ(1));
	t1 = [t94, -t93, 0, 0; t93, t94, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t95 = sin(pkin(6));
	t98 = sin(qJ(1));
	t106 = t98 * t95;
	t97 = sin(qJ(2));
	t105 = t98 * t97;
	t99 = cos(qJ(2));
	t104 = t98 * t99;
	t100 = cos(qJ(1));
	t103 = t100 * t95;
	t102 = t100 * t97;
	t101 = t100 * t99;
	t96 = cos(pkin(6));
	t1 = [-t96 * t105 + t101, -t96 * t104 - t102, t106, t100 * pkin(1) + pkin(8) * t106 + 0; t96 * t102 + t104, t96 * t101 - t105, -t103, t98 * pkin(1) - pkin(8) * t103 + 0; t95 * t97, t95 * t99, t96, t96 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t111 = sin(qJ(2));
	t112 = sin(qJ(1));
	t119 = t112 * t111;
	t113 = cos(qJ(2));
	t118 = t112 * t113;
	t114 = cos(qJ(1));
	t117 = t114 * t111;
	t116 = t114 * t113;
	t115 = pkin(2) * t111 - qJ(3) * t113;
	t110 = cos(pkin(6));
	t109 = sin(pkin(6));
	t108 = t113 * pkin(2) + t111 * qJ(3) + pkin(1);
	t107 = t109 * pkin(8) - t115 * t110;
	t1 = [t112 * t109, t110 * t119 - t116, t110 * t118 + t117, t107 * t112 + t108 * t114 + 0; -t114 * t109, -t110 * t117 - t118, -t110 * t116 + t119, -t107 * t114 + t108 * t112 + 0; t110, -t109 * t111, -t109 * t113, t110 * pkin(8) + t115 * t109 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t123 = sin(pkin(6));
	t128 = cos(qJ(4));
	t140 = t123 * t128;
	t125 = sin(qJ(4));
	t139 = t125 * t123;
	t126 = sin(qJ(2));
	t127 = sin(qJ(1));
	t138 = t127 * t126;
	t129 = cos(qJ(2));
	t137 = t127 * t129;
	t136 = t128 * t129;
	t130 = cos(qJ(1));
	t135 = t130 * t126;
	t134 = t130 * t129;
	t132 = pkin(2) + pkin(9);
	t133 = qJ(3) * t129 - t126 * t132;
	t131 = pkin(3) + pkin(8);
	t124 = cos(pkin(6));
	t122 = t126 * qJ(3) + t132 * t129 + pkin(1);
	t121 = -t124 * t134 + t138;
	t120 = t123 * t131 + t133 * t124;
	t1 = [t127 * t140 + (t124 * t137 + t135) * t125, (t124 * t136 - t139) * t127 + t128 * t135, -t124 * t138 + t134, t120 * t127 + t122 * t130 + 0; t121 * t125 - t130 * t140, t121 * t128 + t130 * t139, t124 * t135 + t137, -t120 * t130 + t122 * t127 + 0; t124 * t128 - t129 * t139, -t123 * t136 - t124 * t125, t123 * t126, -t133 * t123 + t131 * t124 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (65->36), mult. (112->59), div. (0->0), fcn. (146->10), ass. (0->29)
	t149 = sin(qJ(5));
	t151 = sin(qJ(2));
	t168 = t149 * t151;
	t147 = sin(pkin(6));
	t150 = sin(qJ(4));
	t167 = t150 * t147;
	t155 = cos(qJ(2));
	t166 = t150 * t155;
	t153 = cos(qJ(5));
	t165 = t151 * t153;
	t156 = cos(qJ(1));
	t164 = t151 * t156;
	t152 = sin(qJ(1));
	t163 = t152 * t151;
	t154 = cos(qJ(4));
	t162 = t154 * t155;
	t161 = t156 * t155;
	t146 = t150 * pkin(4) - pkin(10) * t154 + qJ(3);
	t157 = pkin(2) + pkin(9);
	t160 = t146 * t155 - t157 * t151;
	t159 = t154 * pkin(4) + t150 * pkin(10) + pkin(3) + pkin(8);
	t148 = cos(pkin(6));
	t143 = t154 * t147 + t148 * t166;
	t158 = t143 * t149 + t148 * t165;
	t145 = -t150 * t168 + t155 * t153;
	t144 = -t147 * t166 + t148 * t154;
	t142 = t146 * t151 + t157 * t155 + pkin(1);
	t141 = t147 * t159 + t160 * t148;
	t1 = [(t143 * t152 + t150 * t164) * t153 + (-t148 * t163 + t161) * t149, t156 * t145 - t158 * t152, (-t148 * t162 + t167) * t152 - t154 * t164, t141 * t152 + t142 * t156 + 0; (-t143 * t153 + t148 * t168) * t156 + t152 * (t155 * t149 + t150 * t165), t152 * t145 + t158 * t156, (t148 * t161 - t163) * t154 - t156 * t167, -t141 * t156 + t142 * t152 + 0; t144 * t153 + t147 * t168, -t144 * t149 + t147 * t165, t147 * t162 + t148 * t150, -t160 * t147 + t159 * t148 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:17:11
	% EndTime: 2020-11-04 22:17:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (91->40), mult. (138->63), div. (0->0), fcn. (172->10), ass. (0->30)
	t179 = sin(qJ(5));
	t181 = sin(qJ(2));
	t197 = t179 * t181;
	t177 = sin(pkin(6));
	t180 = sin(qJ(4));
	t196 = t180 * t177;
	t185 = cos(qJ(2));
	t195 = t180 * t185;
	t183 = cos(qJ(5));
	t194 = t181 * t183;
	t186 = cos(qJ(1));
	t193 = t181 * t186;
	t182 = sin(qJ(1));
	t192 = t182 * t181;
	t184 = cos(qJ(4));
	t191 = t184 * t185;
	t190 = t186 * t185;
	t176 = pkin(5) * t183 + qJ(6) * t179 + pkin(4);
	t171 = -pkin(10) * t184 + t176 * t180 + qJ(3);
	t175 = -t179 * pkin(5) + qJ(6) * t183 - pkin(2) - pkin(9);
	t189 = t171 * t185 + t175 * t181;
	t178 = cos(pkin(6));
	t172 = t184 * t177 + t178 * t195;
	t188 = t172 * t179 + t178 * t194;
	t187 = t180 * pkin(10) + t176 * t184 + pkin(3) + pkin(8);
	t174 = -t180 * t197 + t185 * t183;
	t173 = -t177 * t195 + t178 * t184;
	t170 = t171 * t181 - t175 * t185 + pkin(1);
	t169 = t177 * t187 + t189 * t178;
	t1 = [(t172 * t182 + t180 * t193) * t183 + (-t178 * t192 + t190) * t179, (-t178 * t191 + t196) * t182 - t184 * t193, -t186 * t174 + t188 * t182, t169 * t182 + t170 * t186 + 0; (-t172 * t183 + t178 * t197) * t186 + t182 * (t185 * t179 + t180 * t194), (t178 * t190 - t192) * t184 - t186 * t196, -t182 * t174 - t188 * t186, -t169 * t186 + t170 * t182 + 0; t173 * t183 + t177 * t197, t177 * t191 + t178 * t180, t173 * t179 - t177 * t194, -t189 * t177 + t187 * t178 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end