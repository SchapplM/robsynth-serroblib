% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t97 = cos(qJ(1));
	t96 = sin(qJ(1));
	t1 = [t97, -t96, 0, 0; t96, t97, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t101 = sin(qJ(1));
	t98 = sin(pkin(6));
	t109 = t101 * t98;
	t103 = cos(qJ(1));
	t108 = t103 * t98;
	t100 = sin(qJ(2));
	t107 = t101 * t100;
	t102 = cos(qJ(2));
	t106 = t101 * t102;
	t105 = t103 * t100;
	t104 = t103 * t102;
	t99 = cos(pkin(6));
	t1 = [-t99 * t107 + t104, -t99 * t106 - t105, t109, t103 * pkin(1) + pkin(8) * t109 + 0; t99 * t105 + t106, t99 * t104 - t107, -t108, t101 * pkin(1) - pkin(8) * t108 + 0; t98 * t100, t98 * t102, t99, t99 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t113 = sin(pkin(6));
	t118 = cos(qJ(3));
	t128 = t113 * t118;
	t115 = sin(qJ(3));
	t127 = t115 * t113;
	t116 = sin(qJ(2));
	t126 = t116 * t118;
	t117 = sin(qJ(1));
	t125 = t117 * t116;
	t119 = cos(qJ(2));
	t124 = t117 * t119;
	t120 = cos(qJ(1));
	t123 = t120 * t116;
	t122 = t120 * t119;
	t121 = pkin(2) * t116 - pkin(9) * t119;
	t114 = cos(pkin(6));
	t112 = t119 * pkin(2) + t116 * pkin(9) + pkin(1);
	t111 = t114 * t123 + t124;
	t110 = t113 * pkin(8) - t121 * t114;
	t1 = [(-t114 * t126 + t127) * t117 + t118 * t122, (t114 * t125 - t122) * t115 + t117 * t128, t114 * t124 + t123, t110 * t117 + t112 * t120 + 0; t111 * t118 - t120 * t127, -t111 * t115 - t120 * t128, -t114 * t122 + t125, -t110 * t120 + t112 * t117 + 0; t113 * t126 + t114 * t115, t114 * t118 - t116 * t127, -t113 * t119, t114 * pkin(8) + t121 * t113 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->32), mult. (112->56), div. (0->0), fcn. (146->10), ass. (0->29)
	t139 = sin(pkin(6));
	t145 = cos(qJ(3));
	t156 = t139 * t145;
	t146 = cos(qJ(2));
	t155 = t139 * t146;
	t141 = cos(pkin(6));
	t143 = sin(qJ(2));
	t154 = t141 * t143;
	t153 = t141 * t146;
	t142 = sin(qJ(3));
	t152 = t142 * t139;
	t151 = t143 * t145;
	t150 = t145 * t146;
	t137 = pkin(3) * t145 + qJ(4) * t142 + pkin(2);
	t149 = t146 * pkin(9) - t137 * t143;
	t148 = t142 * pkin(3) - qJ(4) * t145 + pkin(8);
	t147 = cos(qJ(1));
	t144 = sin(qJ(1));
	t140 = cos(pkin(11));
	t138 = sin(pkin(11));
	t136 = t138 * t143 + t140 * t150;
	t135 = t141 * t151 - t152;
	t134 = t139 * t151 + t141 * t142;
	t133 = t138 * t150 - t143 * t140;
	t132 = t143 * pkin(9) + t137 * t146 + pkin(1);
	t131 = t138 * t135 + t140 * t153;
	t130 = -t140 * t135 + t138 * t153;
	t129 = t139 * t148 + t149 * t141;
	t1 = [t130 * t144 + t147 * t136, t131 * t144 - t147 * t133, (-t144 * t154 + t147 * t146) * t142 - t144 * t156, t129 * t144 + t132 * t147 + 0; -t130 * t147 + t144 * t136, -t131 * t147 - t144 * t133, (t144 * t146 + t147 * t154) * t142 + t147 * t156, -t129 * t147 + t132 * t144 + 0; t134 * t140 - t138 * t155, -t134 * t138 - t140 * t155, -t141 * t145 + t143 * t152, -t149 * t139 + t148 * t141 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t170 = sin(pkin(6));
	t176 = cos(qJ(3));
	t190 = t170 * t176;
	t177 = cos(qJ(2));
	t189 = t170 * t177;
	t173 = sin(qJ(3));
	t188 = t173 * t170;
	t174 = sin(qJ(2));
	t187 = t174 * t176;
	t175 = sin(qJ(1));
	t186 = t175 * t174;
	t185 = t175 * t177;
	t178 = cos(qJ(1));
	t184 = t178 * t174;
	t183 = t178 * t177;
	t166 = cos(pkin(11)) * pkin(4) + pkin(3);
	t172 = qJ(4) + pkin(10);
	t159 = t166 * t176 + t172 * t173 + pkin(2);
	t165 = sin(pkin(11)) * pkin(4) + pkin(9);
	t182 = t159 * t174 - t165 * t177;
	t181 = t166 * t173 - t172 * t176 + pkin(8);
	t171 = cos(pkin(6));
	t161 = t171 * t187 - t188;
	t180 = -t161 * t175 + t176 * t183;
	t179 = t161 * t178 + t176 * t185;
	t169 = pkin(11) + qJ(5);
	t168 = cos(t169);
	t167 = sin(t169);
	t163 = t171 * t185 + t184;
	t162 = t171 * t183 - t186;
	t160 = t170 * t187 + t171 * t173;
	t158 = t159 * t177 + t165 * t174 + pkin(1);
	t157 = t181 * t170 - t182 * t171;
	t1 = [t163 * t167 + t180 * t168, t163 * t168 - t180 * t167, (-t171 * t186 + t183) * t173 - t175 * t190, t157 * t175 + t158 * t178 + 0; -t167 * t162 + t179 * t168, -t168 * t162 - t179 * t167, (t171 * t184 + t185) * t173 + t178 * t190, -t157 * t178 + t158 * t175 + 0; t160 * t168 - t167 * t189, -t160 * t167 - t168 * t189, -t171 * t176 + t174 * t188, t182 * t170 + t181 * t171 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:53
	% EndTime: 2020-11-04 22:28:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (128->44), mult. (187->61), div. (0->0), fcn. (241->12), ass. (0->40)
	t210 = sin(pkin(6));
	t216 = cos(qJ(3));
	t230 = t210 * t216;
	t217 = cos(qJ(2));
	t229 = t210 * t217;
	t213 = sin(qJ(3));
	t228 = t213 * t210;
	t214 = sin(qJ(2));
	t227 = t214 * t216;
	t215 = sin(qJ(1));
	t226 = t215 * t214;
	t225 = t215 * t217;
	t218 = cos(qJ(1));
	t224 = t218 * t214;
	t223 = t218 * t217;
	t206 = cos(pkin(11)) * pkin(4) + pkin(3);
	t212 = qJ(4) + pkin(10);
	t199 = t206 * t216 + t212 * t213 + pkin(2);
	t205 = sin(pkin(11)) * pkin(4) + pkin(9);
	t222 = t199 * t214 - t205 * t217;
	t221 = t206 * t213 - t212 * t216 + pkin(8);
	t211 = cos(pkin(6));
	t201 = t211 * t227 - t228;
	t220 = -t201 * t215 + t216 * t223;
	t219 = t201 * t218 + t216 * t225;
	t209 = pkin(11) + qJ(5);
	t208 = cos(t209);
	t207 = sin(t209);
	t203 = t211 * t225 + t224;
	t202 = t211 * t223 - t226;
	t200 = t210 * t227 + t211 * t213;
	t198 = t199 * t217 + t205 * t214 + pkin(1);
	t197 = t200 * t208 - t207 * t229;
	t196 = t200 * t207 + t208 * t229;
	t195 = -t208 * t202 - t219 * t207;
	t194 = t203 * t207 + t220 * t208;
	t193 = -t207 * t202 + t219 * t208;
	t192 = t203 * t208 - t220 * t207;
	t191 = t221 * t210 - t222 * t211;
	t1 = [t194, (-t211 * t226 + t223) * t213 - t215 * t230, -t192, t194 * pkin(5) - t192 * qJ(6) + t191 * t215 + t198 * t218 + 0; t193, (t211 * t224 + t225) * t213 + t218 * t230, -t195, t193 * pkin(5) - t195 * qJ(6) - t191 * t218 + t198 * t215 + 0; t197, -t211 * t216 + t214 * t228, t196, t197 * pkin(5) + t196 * qJ(6) + t222 * t210 + t221 * t211 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end