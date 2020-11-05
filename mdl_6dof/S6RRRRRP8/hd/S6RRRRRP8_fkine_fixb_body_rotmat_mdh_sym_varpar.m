% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t98 = cos(qJ(1));
	t97 = sin(qJ(1));
	t1 = [t98, -t97, 0, 0; t97, t98, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t102 = sin(qJ(1));
	t99 = sin(pkin(6));
	t110 = t102 * t99;
	t104 = cos(qJ(1));
	t109 = t104 * t99;
	t101 = sin(qJ(2));
	t108 = t102 * t101;
	t103 = cos(qJ(2));
	t107 = t102 * t103;
	t106 = t104 * t101;
	t105 = t104 * t103;
	t100 = cos(pkin(6));
	t1 = [-t100 * t108 + t105, -t100 * t107 - t106, t110, t104 * pkin(1) + pkin(8) * t110 + 0; t100 * t106 + t107, t100 * t105 - t108, -t109, t102 * pkin(1) - pkin(8) * t109 + 0; t99 * t101, t99 * t103, t100, t100 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t114 = sin(pkin(6));
	t119 = cos(qJ(3));
	t129 = t114 * t119;
	t116 = sin(qJ(3));
	t128 = t116 * t114;
	t117 = sin(qJ(2));
	t127 = t117 * t119;
	t118 = sin(qJ(1));
	t126 = t118 * t117;
	t120 = cos(qJ(2));
	t125 = t118 * t120;
	t121 = cos(qJ(1));
	t124 = t121 * t117;
	t123 = t121 * t120;
	t122 = pkin(2) * t117 - pkin(9) * t120;
	t115 = cos(pkin(6));
	t113 = t120 * pkin(2) + t117 * pkin(9) + pkin(1);
	t112 = -t115 * t124 - t125;
	t111 = t114 * pkin(8) - t122 * t115;
	t1 = [(-t115 * t127 + t128) * t118 + t119 * t123, (t115 * t126 - t123) * t116 + t118 * t129, t115 * t125 + t124, t111 * t118 + t113 * t121 + 0; -t112 * t119 - t121 * t128, t112 * t116 - t121 * t129, -t115 * t123 + t126, -t111 * t121 + t113 * t118 + 0; t114 * t127 + t115 * t116, t115 * t119 - t117 * t128, -t114 * t120, t115 * pkin(8) + t122 * t114 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t138 = sin(pkin(6));
	t141 = sin(qJ(2));
	t154 = t138 * t141;
	t142 = sin(qJ(1));
	t153 = t138 * t142;
	t144 = cos(qJ(1));
	t152 = t138 * t144;
	t151 = t142 * t141;
	t143 = cos(qJ(2));
	t150 = t142 * t143;
	t149 = t144 * t141;
	t148 = t144 * t143;
	t147 = sin(qJ(3)) * pkin(3) + pkin(8);
	t134 = cos(qJ(3)) * pkin(3) + pkin(2);
	t145 = pkin(10) + pkin(9);
	t146 = t134 * t141 - t145 * t143;
	t139 = cos(pkin(6));
	t137 = qJ(3) + qJ(4);
	t136 = cos(t137);
	t135 = sin(t137);
	t133 = t139 * t149 + t150;
	t132 = t139 * t151 - t148;
	t131 = t134 * t143 + t145 * t141 + pkin(1);
	t130 = t138 * t147 - t146 * t139;
	t1 = [-t132 * t136 + t135 * t153, t132 * t135 + t136 * t153, t139 * t150 + t149, t130 * t142 + t131 * t144 + 0; t133 * t136 - t135 * t152, -t133 * t135 - t136 * t152, -t139 * t148 + t151, -t130 * t144 + t131 * t142 + 0; t139 * t135 + t136 * t154, -t135 * t154 + t139 * t136, -t138 * t143, t146 * t138 + t147 * t139 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (139->56), div. (0->0), fcn. (183->12), ass. (0->36)
	t171 = sin(pkin(6));
	t175 = sin(qJ(2));
	t190 = t171 * t175;
	t176 = sin(qJ(1));
	t189 = t171 * t176;
	t178 = cos(qJ(2));
	t188 = t171 * t178;
	t179 = cos(qJ(1));
	t187 = t171 * t179;
	t186 = t176 * t175;
	t185 = t176 * t178;
	t184 = t179 * t175;
	t183 = t179 * t178;
	t182 = sin(qJ(3)) * pkin(3) + pkin(8);
	t167 = cos(qJ(3)) * pkin(3) + pkin(2);
	t180 = pkin(10) + pkin(9);
	t181 = t167 * t175 - t180 * t178;
	t177 = cos(qJ(5));
	t173 = sin(qJ(5));
	t172 = cos(pkin(6));
	t170 = qJ(3) + qJ(4);
	t169 = cos(t170);
	t168 = sin(t170);
	t166 = t172 * t185 + t184;
	t165 = t172 * t184 + t185;
	t164 = -t172 * t183 + t186;
	t163 = t172 * t186 - t183;
	t162 = t167 * t178 + t180 * t175 + pkin(1);
	t161 = t172 * t168 + t169 * t190;
	t160 = t168 * t190 - t172 * t169;
	t159 = t171 * t182 - t181 * t172;
	t158 = -t163 * t169 + t168 * t189;
	t157 = t165 * t169 - t168 * t187;
	t156 = t165 * t168 + t169 * t187;
	t155 = t163 * t168 + t169 * t189;
	t1 = [t158 * t177 + t166 * t173, -t158 * t173 + t166 * t177, -t155, t158 * pkin(4) - t155 * pkin(11) + t159 * t176 + t162 * t179 + 0; t157 * t177 + t164 * t173, -t157 * t173 + t164 * t177, t156, t157 * pkin(4) + t156 * pkin(11) - t159 * t179 + t162 * t176 + 0; t161 * t177 - t173 * t188, -t161 * t173 - t177 * t188, t160, t161 * pkin(4) + t160 * pkin(11) + t181 * t171 + t182 * t172 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:33
	% EndTime: 2020-11-04 22:44:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (146->43), mult. (198->65), div. (0->0), fcn. (232->14), ass. (0->41)
	t212 = sin(qJ(5));
	t217 = cos(qJ(5));
	t205 = pkin(5) * t217 + qJ(6) * t212 + pkin(4);
	t213 = sin(qJ(4));
	t218 = cos(qJ(4));
	t198 = t213 * pkin(11) + t205 * t218 + pkin(3);
	t214 = sin(qJ(3));
	t194 = t198 * t214 + pkin(8);
	t199 = t218 * pkin(11) - t213 * t205;
	t195 = t199 * t214 + pkin(2);
	t210 = sin(pkin(6));
	t211 = cos(pkin(6));
	t215 = sin(qJ(2));
	t219 = cos(qJ(3));
	t204 = -t212 * pkin(5) + qJ(6) * t217 - pkin(9) - pkin(10);
	t220 = cos(qJ(2));
	t234 = t204 * t220;
	t235 = (t195 * t215 + t234) * t211 + (t211 * t198 * t215 + t210 * t199) * t219 - t194 * t210;
	t233 = t210 * t215;
	t216 = sin(qJ(1));
	t232 = t210 * t216;
	t231 = t210 * t220;
	t221 = cos(qJ(1));
	t230 = t210 * t221;
	t229 = t216 * t215;
	t228 = t216 * t220;
	t227 = t221 * t215;
	t226 = t221 * t220;
	t200 = t211 * t229 - t226;
	t209 = qJ(3) + qJ(4);
	t207 = sin(t209);
	t208 = cos(t209);
	t223 = -t200 * t208 + t207 * t232;
	t202 = t211 * t227 + t228;
	t222 = t202 * t208 - t207 * t230;
	t203 = t211 * t228 + t227;
	t201 = t211 * t226 - t229;
	t197 = t211 * t207 + t208 * t233;
	t193 = t198 * t219 + t195;
	t191 = t193 * t220 - t204 * t215 + pkin(1);
	t1 = [t203 * t212 + t223 * t217, -t200 * t207 - t208 * t232, -t203 * t217 + t223 * t212, t191 * t221 - t235 * t216 + 0; -t212 * t201 + t222 * t217, t202 * t207 + t208 * t230, t217 * t201 + t222 * t212, t191 * t216 + t235 * t221 + 0; t197 * t217 - t212 * t231, t207 * t233 - t211 * t208, t197 * t212 + t217 * t231, pkin(7) + 0 + (-t199 * t219 + t194) * t211 + (t193 * t215 + t234) * t210; 0, 0, 0, 1;];
	Tc_mdh = t1;
end