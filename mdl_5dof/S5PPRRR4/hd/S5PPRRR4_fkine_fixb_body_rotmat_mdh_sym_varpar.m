% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t108 = cos(pkin(10));
	t107 = sin(pkin(10));
	t1 = [t108, -t107, 0, 0; t107, t108, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t110 = sin(pkin(10));
	t111 = sin(pkin(5));
	t118 = t110 * t111;
	t114 = cos(pkin(5));
	t117 = t110 * t114;
	t113 = cos(pkin(10));
	t116 = t113 * t111;
	t115 = t113 * t114;
	t112 = cos(pkin(11));
	t109 = sin(pkin(11));
	t1 = [-t109 * t117 + t113 * t112, -t113 * t109 - t112 * t117, t118, t113 * pkin(1) + qJ(2) * t118 + 0; t109 * t115 + t110 * t112, -t110 * t109 + t112 * t115, -t116, t110 * pkin(1) - qJ(2) * t116 + 0; t111 * t109, t111 * t112, t114, t114 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t130 = sin(pkin(6));
	t147 = pkin(7) * t130;
	t128 = sin(pkin(11));
	t133 = cos(pkin(10));
	t146 = t128 * t133;
	t129 = sin(pkin(10));
	t145 = t129 * t128;
	t135 = cos(pkin(5));
	t144 = t130 * t135;
	t131 = sin(pkin(5));
	t143 = t131 * t130;
	t132 = cos(pkin(11));
	t134 = cos(pkin(6));
	t142 = t132 * t134;
	t141 = t133 * t135;
	t140 = t134 * t131;
	t139 = t135 * t134;
	t125 = -t128 * pkin(2) + t132 * t147;
	t126 = t134 * pkin(7) + qJ(2);
	t138 = t125 * t135 + t131 * t126;
	t137 = cos(qJ(3));
	t136 = sin(qJ(3));
	t124 = t132 * pkin(2) + t128 * t147 + pkin(1);
	t123 = t133 * t132 - t135 * t145;
	t122 = t128 * t141 + t129 * t132;
	t121 = t132 * t139 - t143;
	t120 = -t121 * t129 - t134 * t146;
	t119 = t121 * t133 - t134 * t145;
	t1 = [t120 * t136 + t123 * t137, t120 * t137 - t123 * t136, (t132 * t144 + t140) * t129 + t130 * t146, t124 * t133 + t138 * t129 + 0; t119 * t136 + t122 * t137, t119 * t137 - t122 * t136, (-t132 * t141 + t145) * t130 - t133 * t140, t124 * t129 - t138 * t133 + 0; t136 * t144 + (t128 * t137 + t136 * t142) * t131, t137 * t144 + (-t128 * t136 + t137 * t142) * t131, -t132 * t143 + t139, -t125 * t131 + t126 * t135 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t168 = sin(pkin(11));
	t170 = sin(pkin(6));
	t193 = t168 * t170;
	t174 = cos(pkin(6));
	t192 = t168 * t174;
	t175 = cos(pkin(5));
	t191 = t168 * t175;
	t171 = sin(pkin(5));
	t190 = t171 * t170;
	t172 = cos(pkin(11));
	t189 = t172 * t174;
	t188 = t174 * t171;
	t187 = t175 * t170;
	t186 = t175 * t174;
	t154 = t172 * t186 - t190;
	t169 = sin(pkin(10));
	t173 = cos(pkin(10));
	t150 = t154 * t169 + t173 * t192;
	t157 = -t169 * t191 + t173 * t172;
	t177 = sin(qJ(3));
	t179 = cos(qJ(3));
	t185 = t150 * t177 - t157 * t179;
	t151 = -t154 * t173 + t169 * t192;
	t156 = t169 * t172 + t173 * t191;
	t184 = t151 * t177 - t156 * t179;
	t159 = t172 * t170 * pkin(7) - t168 * pkin(2);
	t165 = t174 * pkin(7) + qJ(2);
	t183 = t159 * t175 + t171 * t165;
	t163 = pkin(3) * t189 + t168 * pkin(8);
	t182 = pkin(3) * t190 - t163 * t175;
	t162 = -t168 * pkin(3) + pkin(8) * t189;
	t181 = pkin(8) * t190 - t162 * t175;
	t180 = t171 * t168 * t179 + (t172 * t188 + t187) * t177;
	t178 = cos(qJ(4));
	t176 = sin(qJ(4));
	t161 = pkin(3) * t192 - t172 * pkin(8);
	t160 = t172 * pkin(3) + pkin(8) * t192;
	t158 = t172 * pkin(2) + pkin(7) * t193 + pkin(1);
	t153 = t172 * t190 - t186;
	t152 = t172 * t187 + t188;
	t149 = t152 * t173 - t169 * t193;
	t148 = t152 * t169 + t173 * t193;
	t1 = [t148 * t176 - t185 * t178, t148 * t178 + t185 * t176, t150 * t179 + t157 * t177, (t173 * t160 - t181 * t169) * t179 + (-t173 * t161 + t182 * t169) * t177 + t183 * t169 + t158 * t173 + 0; -t149 * t176 - t184 * t178, -t178 * t149 + t184 * t176, t151 * t179 + t156 * t177, (t169 * t160 + t181 * t173) * t179 + (-t169 * t161 - t182 * t173) * t177 - t183 * t173 + t158 * t169 + 0; -t176 * t153 + t180 * t178, -t178 * t153 - t180 * t176, -t179 * t187 + (t168 * t177 - t179 * t189) * t171, (-pkin(8) * t187 - t162 * t171) * t179 + (pkin(3) * t187 + t163 * t171) * t177 - t159 * t171 + t165 * t175 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:51
	% EndTime: 2020-11-04 19:55:51
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (167->68), mult. (411->113), div. (0->0), fcn. (536->14), ass. (0->54)
	t223 = sin(pkin(11));
	t225 = sin(pkin(6));
	t250 = t223 * t225;
	t229 = cos(pkin(6));
	t249 = t223 * t229;
	t230 = cos(pkin(5));
	t248 = t223 * t230;
	t226 = sin(pkin(5));
	t247 = t226 * t225;
	t227 = cos(pkin(11));
	t246 = t227 * t229;
	t245 = t229 * t226;
	t244 = t230 * t225;
	t243 = t230 * t229;
	t209 = t227 * t243 - t247;
	t224 = sin(pkin(10));
	t228 = cos(pkin(10));
	t205 = t209 * t224 + t228 * t249;
	t212 = -t224 * t248 + t228 * t227;
	t233 = sin(qJ(3));
	t236 = cos(qJ(3));
	t242 = t205 * t233 - t212 * t236;
	t206 = -t209 * t228 + t224 * t249;
	t211 = t224 * t227 + t228 * t248;
	t241 = t206 * t233 - t211 * t236;
	t214 = t227 * t225 * pkin(7) - t223 * pkin(2);
	t220 = t229 * pkin(7) + qJ(2);
	t240 = t214 * t230 + t226 * t220;
	t218 = pkin(3) * t246 + t223 * pkin(8);
	t239 = pkin(3) * t247 - t218 * t230;
	t217 = -t223 * pkin(3) + pkin(8) * t246;
	t238 = pkin(8) * t247 - t217 * t230;
	t237 = t226 * t223 * t236 + (t227 * t245 + t244) * t233;
	t235 = cos(qJ(4));
	t234 = cos(qJ(5));
	t232 = sin(qJ(4));
	t231 = sin(qJ(5));
	t216 = pkin(3) * t249 - t227 * pkin(8);
	t215 = t227 * pkin(3) + pkin(8) * t249;
	t213 = t227 * pkin(2) + pkin(7) * t250 + pkin(1);
	t208 = t227 * t247 - t243;
	t207 = t227 * t244 + t245;
	t204 = t207 * t228 - t224 * t250;
	t203 = t207 * t224 + t228 * t250;
	t202 = -t236 * t244 + (t223 * t233 - t236 * t246) * t226;
	t201 = -t235 * t208 - t237 * t232;
	t200 = -t232 * t208 + t237 * t235;
	t199 = t206 * t236 + t211 * t233;
	t198 = t205 * t236 + t212 * t233;
	t197 = -t235 * t204 + t241 * t232;
	t196 = t203 * t232 - t242 * t235;
	t195 = -t204 * t232 - t241 * t235;
	t194 = t203 * t235 + t242 * t232;
	t1 = [t196 * t234 + t198 * t231, -t196 * t231 + t198 * t234, -t194, t196 * pkin(4) - t194 * pkin(9) + (t228 * t215 - t238 * t224) * t236 + (-t228 * t216 + t239 * t224) * t233 + t240 * t224 + t213 * t228 + 0; t195 * t234 + t199 * t231, -t195 * t231 + t199 * t234, -t197, t195 * pkin(4) - t197 * pkin(9) + (t224 * t215 + t238 * t228) * t236 + (-t224 * t216 - t239 * t228) * t233 - t240 * t228 + t213 * t224 + 0; t200 * t234 + t202 * t231, -t200 * t231 + t202 * t234, -t201, t200 * pkin(4) - t201 * pkin(9) + (-pkin(8) * t244 - t217 * t226) * t236 + (pkin(3) * t244 + t218 * t226) * t233 - t214 * t226 + t220 * t230 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
end