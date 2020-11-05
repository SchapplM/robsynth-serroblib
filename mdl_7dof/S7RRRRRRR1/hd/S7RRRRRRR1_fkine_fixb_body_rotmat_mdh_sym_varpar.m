% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S7RRRRRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S7RRRRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t116 = cos(qJ(1));
	t115 = sin(qJ(1));
	t1 = [t116, -t115, 0, 0; t115, t116, 0, 0; 0, 0, 1, pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (4->4), div. (0->0), fcn. (12->4), ass. (0->5)
	t120 = cos(qJ(1));
	t119 = cos(qJ(2));
	t118 = sin(qJ(1));
	t117 = sin(qJ(2));
	t1 = [t120 * t119, -t120 * t117, t118, 0; t118 * t119, -t118 * t117, -t120, 0; t117, t119, 0, pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (21->16), div. (0->0), fcn. (34->6), ass. (0->12)
	t122 = sin(qJ(2));
	t123 = sin(qJ(1));
	t131 = t123 * t122;
	t125 = cos(qJ(2));
	t130 = t123 * t125;
	t121 = sin(qJ(3));
	t126 = cos(qJ(1));
	t129 = t126 * t121;
	t128 = t126 * t122;
	t124 = cos(qJ(3));
	t127 = t126 * t124;
	t1 = [-t123 * t121 + t125 * t127, -t123 * t124 - t125 * t129, -t128, -pkin(2) * t128 + 0; t124 * t130 + t129, -t121 * t130 + t127, -t131, -pkin(2) * t131 + 0; t122 * t124, -t122 * t121, t125, t125 * pkin(2) + pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (17->16), mult. (42->30), div. (0->0), fcn. (63->8), ass. (0->17)
	t135 = sin(qJ(2));
	t136 = sin(qJ(1));
	t147 = t135 * t136;
	t138 = cos(qJ(3));
	t146 = t135 * t138;
	t140 = cos(qJ(1));
	t145 = t135 * t140;
	t139 = cos(qJ(2));
	t144 = t136 * t139;
	t137 = cos(qJ(4));
	t143 = t139 * t137;
	t134 = sin(qJ(3));
	t142 = t140 * t134;
	t141 = t140 * t138;
	t133 = sin(qJ(4));
	t132 = -t136 * t134 + t139 * t141;
	t1 = [t132 * t137 + t133 * t145, -t132 * t133 + t137 * t145, -t136 * t138 - t139 * t142, -pkin(2) * t145 + 0; (t135 * t133 + t138 * t143) * t136 + t137 * t142, (-t138 * t144 - t142) * t133 + t137 * t147, -t134 * t144 + t141, -pkin(2) * t147 + 0; -t139 * t133 + t137 * t146, -t133 * t146 - t143, -t135 * t134, t139 * pkin(2) + pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (37->29), mult. (94->57), div. (0->0), fcn. (128->10), ass. (0->28)
	t155 = sin(qJ(3));
	t159 = cos(qJ(4));
	t174 = t155 * t159;
	t161 = cos(qJ(2));
	t173 = t155 * t161;
	t154 = sin(qJ(4));
	t156 = sin(qJ(2));
	t172 = t156 * t154;
	t171 = t156 * t159;
	t157 = sin(qJ(1));
	t170 = t157 * t161;
	t160 = cos(qJ(3));
	t169 = t159 * t160;
	t168 = t161 * t154;
	t167 = t161 * t159;
	t162 = cos(qJ(1));
	t166 = t162 * t155;
	t165 = t162 * t160;
	t164 = pkin(3) * t154 * t155;
	t163 = t160 * t172;
	t158 = cos(qJ(5));
	t153 = sin(qJ(5));
	t152 = t159 * pkin(3) + pkin(2);
	t151 = t160 * t167 + t172;
	t150 = t153 * t174 - t158 * t160;
	t149 = t160 * pkin(3) * t168 - t152 * t156;
	t148 = -t151 * t153 - t158 * t173;
	t1 = [(t151 * t158 - t153 * t173) * t162 - t157 * (t153 * t160 + t158 * t174), t148 * t162 + t157 * t150, (-t157 * t155 + t161 * t165) * t154 - t162 * t171, t149 * t162 - t157 * t164 + 0; (t151 * t157 + t159 * t166) * t158 - t153 * (t155 * t170 - t165), t148 * t157 - t162 * t150, (t160 * t170 + t166) * t154 - t157 * t171, t149 * t157 + t162 * t164 + 0; (-t153 * t155 + t158 * t169) * t156 - t158 * t168, (-t156 * t169 + t168) * t153 - t156 * t155 * t158, t163 + t167, pkin(3) * t163 + t152 * t161 + pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (58->36), mult. (149->69), div. (0->0), fcn. (204->12), ass. (0->34)
	t187 = sin(qJ(5));
	t189 = sin(qJ(3));
	t208 = t187 * t189;
	t188 = sin(qJ(4));
	t207 = t188 * t189;
	t192 = cos(qJ(6));
	t206 = t188 * t192;
	t190 = sin(qJ(2));
	t205 = t190 * t188;
	t193 = cos(qJ(5));
	t204 = t193 * t189;
	t194 = cos(qJ(4));
	t195 = cos(qJ(3));
	t203 = t194 * t195;
	t196 = cos(qJ(2));
	t202 = t196 * t188;
	t201 = t196 * t194;
	t197 = cos(qJ(1));
	t200 = t197 * t207;
	t199 = t195 * t205;
	t198 = t195 * t202;
	t183 = t193 * t203 - t208;
	t177 = t183 * t196 + t193 * t205;
	t191 = sin(qJ(1));
	t186 = sin(qJ(6));
	t185 = t194 * pkin(3) + pkin(2);
	t182 = -t194 * t190 + t198;
	t181 = t187 * t195 + t194 * t204;
	t180 = -t193 * t195 + t194 * t208;
	t179 = pkin(3) * t198 - t185 * t190;
	t178 = (t195 * t201 + t205) * t187 + t196 * t204;
	t176 = t181 * t186 - t189 * t206;
	t175 = -t177 * t186 + t192 * t182;
	t1 = [(t177 * t192 + t186 * t182) * t197 - t191 * (t181 * t192 + t186 * t207), t175 * t197 + t176 * t191, t178 * t197 - t191 * t180, -t191 * pkin(3) * t207 + t179 * t197 + 0; (t177 * t191 + t197 * t181) * t192 + t186 * (t182 * t191 + t200), t175 * t191 - t176 * t197, t178 * t191 + t197 * t180, pkin(3) * t200 + t179 * t191 + 0; (t183 * t190 - t193 * t202) * t192 + (t199 + t201) * t186, (-t183 * t186 + t195 * t206) * t190 + t196 * (t193 * t188 * t186 + t192 * t194), (t190 * t203 - t202) * t187 + t190 * t204, pkin(3) * t199 + t185 * t196 + pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 7
	% From fkine_7_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:50:19
	% EndTime: 2020-11-04 22:50:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (113->51), mult. (285->99), div. (0->0), fcn. (378->14), ass. (0->46)
	t229 = sin(qJ(6));
	t231 = sin(qJ(4));
	t236 = cos(qJ(6));
	t237 = cos(qJ(5));
	t238 = cos(qJ(4));
	t247 = t237 * t238;
	t223 = t229 * t231 + t236 * t247;
	t230 = sin(qJ(5));
	t239 = cos(qJ(3));
	t232 = sin(qJ(3));
	t250 = t232 * t236;
	t218 = -t223 * t239 + t230 * t250;
	t251 = t231 * t237;
	t219 = t229 * t238 - t236 * t251;
	t233 = sin(qJ(2));
	t240 = cos(qJ(2));
	t259 = t218 * t240 + t233 * t219;
	t258 = pkin(4) * t229;
	t227 = t236 * pkin(4) + pkin(3);
	t257 = t227 * t238 + pkin(2);
	t256 = (t230 * t239 + t232 * t247) * t229;
	t246 = t238 * t239;
	t253 = t230 * t232;
	t224 = t237 * t246 - t253;
	t255 = t224 * t240;
	t254 = t230 * t231;
	t252 = t231 * t227;
	t248 = t236 * t239;
	t245 = t239 * t240;
	t244 = t237 * t258;
	t243 = t229 * t251;
	t242 = pkin(4) * t256 - t232 * t252;
	t241 = cos(qJ(1));
	t235 = cos(qJ(7));
	t234 = sin(qJ(1));
	t228 = sin(qJ(7));
	t221 = t237 * t239 - t238 * t253;
	t220 = t230 * t246 + t237 * t232;
	t216 = -t231 * t250 + t256;
	t215 = t220 * t240 + t233 * t254;
	t214 = t232 * t223 + t230 * t248;
	t213 = (-t233 * t251 - t255) * t229 + t236 * (t231 * t245 - t238 * t233);
	t211 = -t214 * t228 + t235 * t221;
	t210 = (t227 * t245 - t233 * t244) * t231 - t255 * t258 - t257 * t233;
	t209 = -t215 * t235 + t259 * t228;
	t1 = [(-t215 * t228 - t235 * t259) * t241 - t234 * (t214 * t235 + t228 * t221), t209 * t241 - t234 * t211, t213 * t241 + t216 * t234, t210 * t241 + t242 * t234 + 0; (t214 * t241 - t234 * t259) * t235 - t228 * (t215 * t234 - t241 * t221), t209 * t234 + t211 * t241, t213 * t234 - t216 * t241, t210 * t234 - t242 * t241 + 0; (-t218 * t233 + t240 * t219) * t235 - (t220 * t233 - t240 * t254) * t228, (t218 * t228 - t235 * t220) * t233 - (t219 * t228 - t235 * t254) * t240, (-t224 * t229 + t231 * t248) * t233 + t240 * (t236 * t238 + t243), (pkin(4) * t243 + t257) * t240 + ((-t238 * t244 + t252) * t239 + t253 * t258) * t233 + 0 + pkin(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
end