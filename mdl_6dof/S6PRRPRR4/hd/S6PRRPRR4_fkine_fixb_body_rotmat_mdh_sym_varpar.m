% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t127 = cos(pkin(11));
	t126 = sin(pkin(11));
	t1 = [t127, -t126, 0, 0; t126, t127, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t128 = sin(pkin(11));
	t129 = sin(pkin(6));
	t137 = t128 * t129;
	t130 = cos(pkin(11));
	t136 = t130 * t129;
	t131 = cos(pkin(6));
	t132 = sin(qJ(2));
	t135 = t131 * t132;
	t133 = cos(qJ(2));
	t134 = t131 * t133;
	t1 = [-t128 * t135 + t130 * t133, -t128 * t134 - t130 * t132, t137, t130 * pkin(1) + pkin(7) * t137 + 0; t128 * t133 + t130 * t135, -t128 * t132 + t130 * t134, -t136, t128 * pkin(1) - pkin(7) * t136 + 0; t129 * t132, t129 * t133, t131, t131 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:28
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t141 = sin(pkin(6));
	t154 = t141 * pkin(7);
	t140 = sin(pkin(11));
	t143 = cos(pkin(6));
	t153 = t140 * t143;
	t144 = sin(qJ(3));
	t152 = t141 * t144;
	t146 = cos(qJ(3));
	t151 = t141 * t146;
	t142 = cos(pkin(11));
	t150 = t142 * t143;
	t145 = sin(qJ(2));
	t149 = t143 * t145;
	t147 = cos(qJ(2));
	t148 = t143 * t147;
	t139 = t140 * t147 + t142 * t149;
	t138 = t140 * t149 - t142 * t147;
	t1 = [-t138 * t146 + t140 * t152, t138 * t144 + t140 * t151, t140 * t148 + t142 * t145, (t142 * pkin(2) + pkin(8) * t153) * t147 + (-pkin(2) * t153 + t142 * pkin(8)) * t145 + t140 * t154 + t142 * pkin(1) + 0; t139 * t146 - t142 * t152, -t139 * t144 - t142 * t151, t140 * t145 - t142 * t148, (t140 * pkin(2) - pkin(8) * t150) * t147 + (pkin(2) * t150 + t140 * pkin(8)) * t145 - t142 * t154 + t140 * pkin(1) + 0; t143 * t144 + t145 * t151, t143 * t146 - t145 * t152, -t141 * t147, t143 * pkin(7) + qJ(1) + 0 + (pkin(2) * t145 - pkin(8) * t147) * t141; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t158 = sin(pkin(6));
	t160 = cos(pkin(6));
	t162 = sin(qJ(2));
	t164 = cos(qJ(2));
	t161 = sin(qJ(3));
	t163 = cos(qJ(3));
	t168 = pkin(3) * t163 + qJ(4) * t161 + pkin(2);
	t166 = -pkin(8) * t164 + t168 * t162;
	t167 = t161 * pkin(3) - qJ(4) * t163 + pkin(7);
	t176 = t167 * t158 - t166 * t160;
	t172 = t160 * t162;
	t171 = t160 * t164;
	t170 = t161 * t158;
	t169 = t163 * t158;
	t165 = pkin(8) * t162 + t168 * t164 + pkin(1);
	t159 = cos(pkin(11));
	t157 = sin(pkin(11));
	t156 = t157 * t164 + t159 * t172;
	t155 = t157 * t172 - t159 * t164;
	t1 = [-t155 * t163 + t157 * t170, t157 * t171 + t159 * t162, -t155 * t161 - t157 * t169, t176 * t157 + t165 * t159 + 0; t156 * t163 - t159 * t170, t157 * t162 - t159 * t171, t156 * t161 + t159 * t169, t165 * t157 - t176 * t159 + 0; t160 * t161 + t162 * t169, -t158 * t164, -t160 * t163 + t162 * t170, t166 * t158 + t167 * t160 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:28
	% EndTime: 2020-11-04 21:11:29
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (73->33), mult. (144->54), div. (0->0), fcn. (181->10), ass. (0->31)
	t185 = sin(pkin(6));
	t187 = cos(pkin(6));
	t190 = sin(qJ(2));
	t193 = cos(qJ(2));
	t194 = pkin(8) - pkin(9);
	t189 = sin(qJ(3));
	t192 = cos(qJ(3));
	t195 = pkin(3) + pkin(4);
	t200 = qJ(4) * t189 + t195 * t192 + pkin(2);
	t197 = t200 * t190 - t194 * t193;
	t199 = qJ(4) * t192 - t195 * t189 - pkin(7);
	t209 = t199 * t185 + t197 * t187;
	t188 = sin(qJ(5));
	t206 = t185 * t188;
	t205 = t185 * t190;
	t191 = cos(qJ(5));
	t204 = t185 * t191;
	t203 = t187 * t190;
	t202 = t187 * t193;
	t184 = sin(pkin(11));
	t186 = cos(pkin(11));
	t180 = t184 * t203 - t186 * t193;
	t198 = t180 * t191 + t184 * t206;
	t196 = t194 * t190 + t200 * t193 + pkin(1);
	t183 = t187 * t189 + t192 * t205;
	t182 = -t187 * t192 + t189 * t205;
	t181 = t184 * t193 + t186 * t203;
	t179 = t181 * t191 + t186 * t206;
	t178 = -t188 * t181 + t186 * t204;
	t177 = -t188 * t180 + t184 * t204;
	t1 = [t177 * t189 - t198 * t192, -t177 * t192 - t198 * t189, -t184 * t202 - t186 * t190, -t209 * t184 + t196 * t186 + 0; -t178 * t189 + t179 * t192, t178 * t192 + t179 * t189, -t184 * t190 + t186 * t202, t196 * t184 + t209 * t186 + 0; t182 * t188 + t183 * t191, t182 * t191 - t183 * t188, t185 * t193, t197 * t185 - t199 * t187 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:29
	% EndTime: 2020-11-04 21:11:29
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (126->52), mult. (225->85), div. (0->0), fcn. (283->12), ass. (0->41)
	t229 = sin(qJ(5));
	t233 = cos(qJ(5));
	t221 = t233 * pkin(5) + t229 * pkin(10) + pkin(3) + pkin(4);
	t222 = -t229 * pkin(5) + t233 * pkin(10) - qJ(4);
	t225 = sin(pkin(6));
	t227 = cos(pkin(6));
	t231 = sin(qJ(2));
	t234 = cos(qJ(3));
	t230 = sin(qJ(3));
	t240 = t221 * t230 + pkin(7);
	t241 = t222 * t230 - pkin(2);
	t235 = cos(qJ(2));
	t236 = pkin(8) - pkin(9);
	t242 = t236 * t235;
	t248 = t227 * t231;
	t256 = (t241 * t231 + t242) * t227 + t240 * t225 - (t221 * t248 - t225 * t222) * t234;
	t224 = sin(pkin(11));
	t252 = t224 * t235;
	t251 = t225 * t231;
	t250 = t225 * t235;
	t226 = cos(pkin(11));
	t249 = t226 * t235;
	t247 = t227 * t235;
	t246 = t229 * t225;
	t245 = t229 * t231;
	t244 = t233 * t225;
	t243 = t233 * t227;
	t239 = t230 * (t225 * t245 + t243) + (-t229 * t227 + t231 * t244) * t234;
	t238 = t221 * t234 - t241;
	t237 = t236 * t231 + t238 * t235 + pkin(1);
	t232 = cos(qJ(6));
	t228 = sin(qJ(6));
	t220 = t230 * t229 + t234 * t233;
	t217 = t224 * t247 + t226 * t231;
	t216 = t226 * t248 + t252;
	t215 = t224 * t231 - t226 * t247;
	t214 = t224 * t248 - t249;
	t212 = (t231 * t243 + t246) * t234 - t230 * (-t227 * t245 + t244);
	t211 = -t224 * t212 + t220 * t249;
	t210 = t226 * t212 + t220 * t252;
	t1 = [t211 * t232 - t228 * t217, -t211 * t228 - t232 * t217, (-t229 * t214 + t224 * t244) * t234 + (t214 * t233 + t224 * t246) * t230, t256 * t224 + t237 * t226 + 0; t210 * t232 - t228 * t215, -t210 * t228 - t232 * t215, (t229 * t216 - t226 * t244) * t234 - (t216 * t233 + t226 * t246) * t230, t237 * t224 - t256 * t226 + 0; t228 * t250 + t239 * t232, -t239 * t228 + t232 * t250, (t227 * t230 + t234 * t251) * t229 - (-t227 * t234 + t230 * t251) * t233, qJ(1) + 0 + (t238 * t231 - t242) * t225 + (t222 * t234 + t240) * t227; 0, 0, 0, 1;];
	Tc_mdh = t1;
end