% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:59
	% EndTime: 2020-11-04 22:32:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:59
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t116 = cos(qJ(1));
	t115 = sin(qJ(1));
	t1 = [t116, -t115, 0, 0; t115, t116, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:00
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t117 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t120 * t117;
	t119 = sin(qJ(2));
	t127 = t120 * t119;
	t121 = cos(qJ(2));
	t126 = t120 * t121;
	t122 = cos(qJ(1));
	t125 = t122 * t117;
	t124 = t122 * t119;
	t123 = t122 * t121;
	t118 = cos(pkin(6));
	t1 = [-t118 * t127 + t123, -t118 * t126 - t124, t128, t122 * pkin(1) + pkin(8) * t128 + 0; t118 * t124 + t126, t118 * t123 - t127, -t125, t120 * pkin(1) - pkin(8) * t125 + 0; t117 * t119, t117 * t121, t118, t118 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:00
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t132 = sin(pkin(6));
	t137 = cos(qJ(3));
	t147 = t132 * t137;
	t134 = sin(qJ(3));
	t146 = t134 * t132;
	t135 = sin(qJ(2));
	t145 = t135 * t137;
	t136 = sin(qJ(1));
	t144 = t136 * t135;
	t138 = cos(qJ(2));
	t143 = t136 * t138;
	t139 = cos(qJ(1));
	t142 = t139 * t135;
	t141 = t139 * t138;
	t140 = pkin(2) * t135 - pkin(9) * t138;
	t133 = cos(pkin(6));
	t131 = t138 * pkin(2) + t135 * pkin(9) + pkin(1);
	t130 = -t133 * t142 - t143;
	t129 = t132 * pkin(8) - t140 * t133;
	t1 = [(-t133 * t145 + t146) * t136 + t137 * t141, (t133 * t144 - t141) * t134 + t136 * t147, t133 * t143 + t142, t129 * t136 + t131 * t139 + 0; -t130 * t137 - t139 * t146, t130 * t134 - t139 * t147, -t133 * t141 + t144, -t129 * t139 + t131 * t136 + 0; t132 * t145 + t133 * t134, t133 * t137 - t135 * t146, -t132 * t138, t133 * pkin(8) + t140 * t132 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:00
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t152 = sin(pkin(6));
	t157 = cos(qJ(3));
	t168 = t152 * t157;
	t154 = sin(qJ(3));
	t167 = t154 * t152;
	t155 = sin(qJ(2));
	t166 = t155 * t157;
	t156 = sin(qJ(1));
	t165 = t156 * t155;
	t158 = cos(qJ(2));
	t164 = t156 * t158;
	t159 = cos(qJ(1));
	t163 = t159 * t155;
	t162 = t159 * t158;
	t151 = t157 * pkin(3) + qJ(4) * t154 + pkin(2);
	t161 = t158 * pkin(9) - t151 * t155;
	t160 = t154 * pkin(3) - qJ(4) * t157 + pkin(8);
	t153 = cos(pkin(6));
	t150 = t153 * t163 + t164;
	t149 = t155 * pkin(9) + t151 * t158 + pkin(1);
	t148 = t152 * t160 + t161 * t153;
	t1 = [(-t153 * t166 + t167) * t156 + t157 * t162, t153 * t164 + t163, (-t153 * t165 + t162) * t154 - t156 * t168, t148 * t156 + t149 * t159 + 0; t150 * t157 - t159 * t167, -t153 * t162 + t165, t150 * t154 + t159 * t168, -t148 * t159 + t149 * t156 + 0; t152 * t166 + t153 * t154, -t152 * t158, -t153 * t157 + t155 * t167, -t161 * t152 + t160 * t153 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:00
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->33), mult. (116->53), div. (0->0), fcn. (153->10), ass. (0->30)
	t179 = sin(pkin(6));
	t183 = sin(qJ(2));
	t198 = t179 * t183;
	t180 = cos(pkin(6));
	t197 = t180 * t183;
	t184 = sin(qJ(1));
	t187 = cos(qJ(2));
	t196 = t184 * t187;
	t188 = cos(qJ(1));
	t195 = t187 * t188;
	t181 = sin(qJ(5));
	t185 = cos(qJ(5));
	t171 = -t179 * t185 + t181 * t197;
	t174 = t179 * t181 + t185 * t197;
	t182 = sin(qJ(3));
	t186 = cos(qJ(3));
	t194 = t171 * t186 - t182 * t174;
	t193 = t182 * t171 + t174 * t186;
	t190 = pkin(3) + pkin(4);
	t177 = qJ(4) * t182 + t190 * t186 + pkin(2);
	t189 = pkin(9) - pkin(10);
	t192 = t177 * t183 - t189 * t187;
	t191 = qJ(4) * t186 - t190 * t182 - pkin(8);
	t176 = t182 * t181 + t186 * t185;
	t175 = t186 * t181 - t185 * t182;
	t173 = t180 * t182 + t186 * t198;
	t172 = -t180 * t186 + t182 * t198;
	t170 = t177 * t187 + t189 * t183 + pkin(1);
	t169 = -t191 * t179 - t192 * t180;
	t1 = [t176 * t195 - t193 * t184, -t175 * t195 + t194 * t184, -t180 * t196 - t188 * t183, t169 * t184 + t170 * t188 + 0; t176 * t196 + t193 * t188, -t175 * t196 - t194 * t188, t180 * t195 - t184 * t183, -t169 * t188 + t170 * t184 + 0; t172 * t181 + t173 * t185, t172 * t185 - t173 * t181, t179 * t187, t192 * t179 - t191 * t180 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:00
	% EndTime: 2020-11-04 22:33:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (126->48), mult. (208->79), div. (0->0), fcn. (266->12), ass. (0->39)
	t221 = sin(qJ(5));
	t226 = cos(qJ(5));
	t214 = pkin(5) * t226 + pkin(11) * t221 + pkin(3) + pkin(4);
	t227 = cos(qJ(3));
	t215 = -t221 * pkin(5) + pkin(11) * t226 - qJ(4);
	t222 = sin(qJ(3));
	t246 = -t215 * t222 + pkin(2);
	t250 = t214 * t227 + t246;
	t218 = sin(pkin(6));
	t219 = cos(pkin(6));
	t223 = sin(qJ(2));
	t242 = t221 * t223;
	t203 = -t218 * t226 + t219 * t242;
	t241 = t223 * t226;
	t206 = t218 * t221 + t219 * t241;
	t199 = t222 * t203 + t206 * t227;
	t213 = t222 * t221 + t227 * t226;
	t224 = sin(qJ(1));
	t228 = cos(qJ(2));
	t229 = cos(qJ(1));
	t239 = t228 * t229;
	t249 = t199 * t224 - t213 * t239;
	t240 = t224 * t228;
	t248 = t199 * t229 + t213 * t240;
	t207 = t214 * t222 + pkin(8);
	t230 = pkin(9) - pkin(10);
	t238 = t228 * t230;
	t247 = (-t223 * t246 + t238) * t219 - (t219 * t214 * t223 - t218 * t215) * t227 + t207 * t218;
	t244 = t218 * t223;
	t243 = t218 * t228;
	t234 = t203 * t227 - t222 * t206;
	t233 = t222 * (t218 * t242 + t219 * t226) + (t218 * t241 - t219 * t221) * t227;
	t231 = t230 * t223 + t250 * t228 + pkin(1);
	t225 = cos(qJ(6));
	t220 = sin(qJ(6));
	t212 = t227 * t221 - t226 * t222;
	t209 = t219 * t240 + t229 * t223;
	t208 = t219 * t239 - t224 * t223;
	t1 = [-t220 * t209 - t249 * t225, -t225 * t209 + t249 * t220, t212 * t239 - t234 * t224, t247 * t224 + t231 * t229 + 0; t208 * t220 + t248 * t225, t208 * t225 - t248 * t220, t212 * t240 + t234 * t229, t231 * t224 - t247 * t229 + 0; t220 * t243 + t233 * t225, -t233 * t220 + t225 * t243, (t219 * t222 + t227 * t244) * t221 - (-t219 * t227 + t222 * t244) * t226, pkin(7) + 0 + (t250 * t223 - t238) * t218 + (t215 * t227 + t207) * t219; 0, 0, 0, 1;];
	Tc_mdh = t1;
end