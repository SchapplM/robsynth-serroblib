% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:03
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t113 = cos(qJ(1));
	t112 = sin(qJ(1));
	t1 = [t113, -t112, 0, 0; t112, t113, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t114 = sin(pkin(6));
	t117 = sin(qJ(1));
	t125 = t117 * t114;
	t116 = sin(qJ(2));
	t124 = t117 * t116;
	t118 = cos(qJ(2));
	t123 = t117 * t118;
	t119 = cos(qJ(1));
	t122 = t119 * t114;
	t121 = t119 * t116;
	t120 = t119 * t118;
	t115 = cos(pkin(6));
	t1 = [-t115 * t124 + t120, -t115 * t123 - t121, t125, t119 * pkin(1) + pkin(8) * t125 + 0; t115 * t121 + t123, t115 * t120 - t124, -t122, t117 * pkin(1) - pkin(8) * t122 + 0; t114 * t116, t114 * t118, t115, t115 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t145 = pkin(2) * sin(qJ(2));
	t138 = qJ(2) + pkin(11);
	t144 = cos(qJ(1));
	t143 = sin(qJ(1));
	t141 = pkin(8) + qJ(3);
	t140 = cos(pkin(6));
	t139 = sin(pkin(6));
	t137 = cos(t138);
	t136 = sin(t138);
	t135 = pkin(6) - t138;
	t134 = pkin(6) + t138;
	t133 = cos(qJ(2)) * pkin(2) + pkin(1);
	t132 = cos(t135);
	t131 = cos(t134);
	t130 = sin(t135);
	t129 = sin(t134);
	t128 = t131 + t132;
	t127 = -t129 + t130;
	t126 = -t139 * t141 + t140 * t145;
	t1 = [t144 * t137 + t143 * t127 / 0.2e1, -t144 * t136 - t143 * t128 / 0.2e1, t143 * t139, -t126 * t143 + t144 * t133 + 0; t143 * t137 - t144 * t127 / 0.2e1, -t143 * t136 + t144 * t128 / 0.2e1, -t144 * t139, t144 * t126 + t143 * t133 + 0; t132 / 0.2e1 - t131 / 0.2e1, t130 / 0.2e1 + t129 / 0.2e1, t140, t139 * t145 + t140 * t141 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (73->27), mult. (61->34), div. (0->0), fcn. (72->14), ass. (0->26)
	t160 = qJ(2) + pkin(11);
	t161 = sin(pkin(11));
	t163 = cos(pkin(11));
	t150 = t163 * pkin(3) + qJ(4) * t161 + pkin(2);
	t151 = -t161 * pkin(3) + qJ(4) * t163;
	t166 = sin(qJ(2));
	t168 = cos(qJ(2));
	t170 = t150 * t166 - t151 * t168;
	t169 = cos(qJ(1));
	t167 = sin(qJ(1));
	t165 = pkin(8) + qJ(3);
	t164 = cos(pkin(6));
	t162 = sin(pkin(6));
	t159 = cos(t160);
	t158 = sin(t160);
	t157 = pkin(6) - t160;
	t156 = pkin(6) + t160;
	t155 = cos(t157);
	t154 = cos(t156);
	t153 = sin(t157);
	t152 = sin(t156);
	t149 = t154 + t155;
	t148 = -t152 + t153;
	t147 = t150 * t168 + t151 * t166 + pkin(1);
	t146 = t162 * t165 - t170 * t164;
	t1 = [t167 * t162, -t169 * t159 - t167 * t148 / 0.2e1, t169 * t158 + t167 * t149 / 0.2e1, t146 * t167 + t147 * t169 + 0; -t169 * t162, -t167 * t159 + t169 * t148 / 0.2e1, t167 * t158 - t169 * t149 / 0.2e1, -t146 * t169 + t147 * t167 + 0; t164, -t155 / 0.2e1 + t154 / 0.2e1, -t153 / 0.2e1 - t152 / 0.2e1, t170 * t162 + t164 * t165 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (98->36), mult. (127->54), div. (0->0), fcn. (115->16), ass. (0->35)
	t194 = sin(qJ(1));
	t206 = t194 / 0.2e1;
	t189 = sin(pkin(6));
	t192 = sin(qJ(5));
	t205 = t189 * t192;
	t195 = cos(qJ(5));
	t204 = t189 * t195;
	t197 = cos(qJ(1));
	t203 = t189 * t197;
	t187 = qJ(2) + pkin(11);
	t185 = cos(t187);
	t202 = t194 * t185;
	t201 = t197 * t185;
	t182 = pkin(6) + t187;
	t178 = sin(t182);
	t183 = pkin(6) - t187;
	t179 = sin(t183);
	t200 = t179 / 0.2e1 - t178 / 0.2e1;
	t199 = cos(t183) / 0.2e1 - cos(t182) / 0.2e1;
	t188 = sin(pkin(11));
	t190 = cos(pkin(11));
	t176 = t190 * pkin(3) + qJ(4) * t188 + pkin(2);
	t177 = -t188 * pkin(3) + qJ(4) * t190;
	t193 = sin(qJ(2));
	t196 = cos(qJ(2));
	t198 = t176 * t193 - t177 * t196;
	t191 = cos(pkin(6));
	t186 = qJ(3) + pkin(4) + pkin(8);
	t184 = sin(t187);
	t175 = -t178 + t179;
	t174 = t197 * t184 + t191 * t202;
	t173 = t194 * t184 - t191 * t201;
	t172 = 0.2e1 * t176 * t196 + 0.2e1 * t177 * t193 + (2 * pkin(1));
	t171 = t189 * t186 - t198 * t191;
	t1 = [t174 * t192 + t194 * t204, t174 * t195 - t194 * t205, t175 * t206 + t201, t171 * t194 + t172 * t197 / 0.2e1 + 0 + (t200 * t194 + t201) * pkin(9); t173 * t192 - t195 * t203, t173 * t195 + t192 * t203, t202 - t197 * t175 / 0.2e1, -t171 * t197 + t172 * t206 + 0 + (-t200 * t197 + t202) * pkin(9); -t185 * t205 + t191 * t195, -t185 * t204 - t191 * t192, t199, t199 * pkin(9) + t186 * t191 + t198 * t189 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:16
	% EndTime: 2020-11-04 22:03:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (118->55), mult. (174->81), div. (0->0), fcn. (202->16), ass. (0->41)
	t217 = sin(pkin(11));
	t247 = qJ(4) * t217 + pkin(2);
	t216 = qJ(2) + pkin(11);
	t212 = sin(t216);
	t218 = sin(pkin(6));
	t246 = t212 * t218;
	t213 = cos(t216);
	t220 = cos(pkin(6));
	t245 = t213 * t220;
	t222 = sin(qJ(5));
	t244 = t218 * t222;
	t226 = cos(qJ(5));
	t243 = t218 * t226;
	t221 = sin(qJ(6));
	t228 = cos(qJ(1));
	t242 = t221 * t228;
	t241 = t222 * t228;
	t224 = sin(qJ(1));
	t240 = t224 * t221;
	t225 = cos(qJ(6));
	t239 = t225 * t224;
	t238 = t228 * t225;
	t237 = t222 * t240;
	t236 = t222 * t239;
	t235 = t224 * t243;
	t234 = t228 * t243;
	t233 = t221 * t241;
	t232 = t222 * t238;
	t231 = t222 * pkin(5) - t226 * pkin(10);
	t230 = pkin(5) * t226 + t222 * pkin(10) + pkin(4) + pkin(8) + qJ(3);
	t229 = pkin(3) + pkin(9);
	t227 = cos(qJ(2));
	t223 = sin(qJ(2));
	t219 = cos(pkin(11));
	t215 = qJ(4) * t219;
	t211 = -t213 * t244 + t220 * t226;
	t210 = -t217 * t229 + t231 * t219 + t215;
	t209 = t231 * t217 + t229 * t219 + t247;
	t208 = t209 * t227 + t210 * t223 + pkin(1);
	t207 = t218 * t230 + (-t209 * t223 + t210 * t227) * t220;
	t1 = [(t220 * t236 + t242) * t213 + (-t220 * t240 + t232) * t212 + t225 * t235, (-t220 * t237 + t238) * t213 + (-t220 * t239 - t233) * t212 - t221 * t235, t224 * t244 - (t228 * t212 + t224 * t245) * t226, t207 * t224 + t208 * t228 + 0; (-t220 * t232 + t240) * t213 + (t220 * t242 + t236) * t212 - t225 * t234, (t220 * t233 + t239) * t213 + (t220 * t238 - t237) * t212 + t221 * t234, -t218 * t241 + (-t224 * t212 + t228 * t245) * t226, -t207 * t228 + t208 * t224 + 0; t211 * t225 + t221 * t246, -t211 * t221 + t225 * t246, t213 * t243 + t220 * t222, pkin(7) + 0 + (cos(pkin(6) - t216) / 0.2e1 - cos(pkin(6) + t216) / 0.2e1) * pkin(9) + t230 * t220 + (-t231 * t213 + (t219 * pkin(3) + t247) * t223 - (-t217 * pkin(3) + t215) * t227) * t218; 0, 0, 0, 1;];
	Tc_mdh = t1;
end