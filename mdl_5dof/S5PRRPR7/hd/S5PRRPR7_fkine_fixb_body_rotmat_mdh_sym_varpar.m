% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t104 = cos(pkin(9));
	t103 = sin(pkin(9));
	t1 = [t104, -t103, 0, 0; t103, t104, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t105 = sin(pkin(9));
	t106 = sin(pkin(5));
	t114 = t105 * t106;
	t107 = cos(pkin(9));
	t113 = t107 * t106;
	t108 = cos(pkin(5));
	t109 = sin(qJ(2));
	t112 = t108 * t109;
	t110 = cos(qJ(2));
	t111 = t108 * t110;
	t1 = [-t105 * t112 + t107 * t110, -t105 * t111 - t107 * t109, t114, t107 * pkin(1) + pkin(6) * t114 + 0; t105 * t110 + t107 * t112, -t105 * t109 + t107 * t111, -t113, t105 * pkin(1) - pkin(6) * t113 + 0; t106 * t109, t106 * t110, t108, t108 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t118 = sin(pkin(5));
	t131 = t118 * pkin(6);
	t117 = sin(pkin(9));
	t120 = cos(pkin(5));
	t130 = t117 * t120;
	t121 = sin(qJ(3));
	t129 = t118 * t121;
	t123 = cos(qJ(3));
	t128 = t118 * t123;
	t119 = cos(pkin(9));
	t127 = t119 * t120;
	t122 = sin(qJ(2));
	t126 = t120 * t122;
	t124 = cos(qJ(2));
	t125 = t120 * t124;
	t116 = t117 * t124 + t119 * t126;
	t115 = t117 * t126 - t119 * t124;
	t1 = [-t115 * t123 + t117 * t129, t115 * t121 + t117 * t128, t117 * t125 + t119 * t122, (t119 * pkin(2) + pkin(7) * t130) * t124 + (-pkin(2) * t130 + t119 * pkin(7)) * t122 + t117 * t131 + t119 * pkin(1) + 0; t116 * t123 - t119 * t129, -t116 * t121 - t119 * t128, t117 * t122 - t119 * t125, (t117 * pkin(2) - pkin(7) * t127) * t124 + (pkin(2) * t127 + t117 * pkin(7)) * t122 - t119 * t131 + t117 * pkin(1) + 0; t120 * t121 + t122 * t128, t120 * t123 - t122 * t129, -t118 * t124, t120 * pkin(6) + qJ(1) + 0 + (pkin(2) * t122 - pkin(7) * t124) * t118; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->39), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->29)
	t144 = sin(pkin(5));
	t159 = t144 * pkin(6);
	t143 = sin(pkin(9));
	t147 = cos(pkin(5));
	t158 = t143 * t147;
	t148 = sin(qJ(3));
	t157 = t144 * t148;
	t150 = cos(qJ(3));
	t156 = t144 * t150;
	t151 = cos(qJ(2));
	t155 = t144 * t151;
	t146 = cos(pkin(9));
	t154 = t146 * t147;
	t149 = sin(qJ(2));
	t153 = t147 * t149;
	t152 = t147 * t151;
	t145 = cos(pkin(10));
	t142 = sin(pkin(10));
	t141 = t147 * t148 + t149 * t156;
	t140 = -t147 * t150 + t149 * t157;
	t139 = t143 * t152 + t146 * t149;
	t138 = t143 * t151 + t146 * t153;
	t137 = t143 * t149 - t146 * t152;
	t136 = t143 * t153 - t146 * t151;
	t135 = -t136 * t150 + t143 * t157;
	t134 = t138 * t150 - t146 * t157;
	t133 = t138 * t148 + t146 * t156;
	t132 = t136 * t148 + t143 * t156;
	t1 = [t135 * t145 + t139 * t142, -t135 * t142 + t139 * t145, -t132, t135 * pkin(3) - t132 * qJ(4) + (t146 * pkin(2) + pkin(7) * t158) * t151 + (-pkin(2) * t158 + t146 * pkin(7)) * t149 + t143 * t159 + t146 * pkin(1) + 0; t134 * t145 + t137 * t142, -t134 * t142 + t137 * t145, t133, t134 * pkin(3) + t133 * qJ(4) + (t143 * pkin(2) - pkin(7) * t154) * t151 + (pkin(2) * t154 + t143 * pkin(7)) * t149 - t146 * t159 + t143 * pkin(1) + 0; t141 * t145 - t142 * t155, -t141 * t142 - t145 * t155, t140, t141 * pkin(3) + t147 * pkin(6) + t140 * qJ(4) + qJ(1) + 0 + (pkin(2) * t149 - pkin(7) * t151) * t144; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:45
	% EndTime: 2020-11-04 20:04:45
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (104->47), mult. (235->81), div. (0->0), fcn. (290->12), ass. (0->42)
	t172 = sin(pkin(5));
	t175 = cos(pkin(5));
	t170 = sin(pkin(10));
	t173 = cos(pkin(10));
	t168 = -t170 * pkin(4) + t173 * pkin(8) - pkin(7);
	t178 = sin(qJ(2));
	t181 = cos(qJ(2));
	t169 = t173 * pkin(4) + t170 * pkin(8) + pkin(3);
	t177 = sin(qJ(3));
	t180 = cos(qJ(3));
	t190 = qJ(4) * t177 + t169 * t180 + pkin(2);
	t183 = t168 * t181 + t190 * t178;
	t189 = qJ(4) * t180 - t169 * t177 - pkin(6);
	t210 = t189 * t172 + t183 * t175;
	t171 = sin(pkin(9));
	t174 = cos(pkin(9));
	t203 = t170 * t175;
	t193 = t171 * t203;
	t199 = t173 * t180;
	t202 = t170 * t178;
	t195 = t178 * t180;
	t200 = t172 * t177;
	t206 = (t175 * t195 - t200) * t173;
	t209 = t171 * t206 - (t174 * t199 + t193) * t181 - t174 * t202;
	t192 = t174 * t203;
	t208 = t174 * t206 + (t171 * t199 - t192) * t181 + t171 * t202;
	t201 = t170 * t181;
	t198 = t175 * t177;
	t197 = t177 * t178;
	t196 = t177 * t181;
	t194 = t181 * t173;
	t191 = t180 * t201;
	t182 = -t168 * t178 + t190 * t181 + pkin(1);
	t179 = cos(qJ(5));
	t176 = sin(qJ(5));
	t167 = -t172 * t197 + t175 * t180;
	t166 = t172 * t180 + t175 * t197;
	t165 = t170 * t200 - t175 * t194;
	t162 = t171 * t166 - t174 * t196;
	t161 = t174 * t166 + t171 * t196;
	t160 = t173 * t198 + (t173 * t195 - t201) * t172;
	t1 = [-t176 * t162 - t209 * t179, -t179 * t162 + t209 * t176, (-t173 * t174 - t180 * t193) * t178 + t174 * t191 + t171 * t165, -t171 * t210 + t182 * t174 + 0; t176 * t161 + t208 * t179, t161 * t179 - t208 * t176, (-t173 * t171 + t180 * t192) * t178 + t171 * t191 - t174 * t165, t182 * t171 + t174 * t210 + 0; t160 * t179 - t176 * t167, -t160 * t176 - t179 * t167, (t172 * t195 + t198) * t170 + t172 * t194, t183 * t172 - t189 * t175 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end