% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t100 = cos(qJ(1));
	t99 = sin(qJ(1));
	t1 = [t100, -t99, 0, 0; t99, t100, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t101 = sin(pkin(6));
	t104 = sin(qJ(1));
	t112 = t104 * t101;
	t103 = sin(qJ(2));
	t111 = t104 * t103;
	t105 = cos(qJ(2));
	t110 = t104 * t105;
	t106 = cos(qJ(1));
	t109 = t106 * t101;
	t108 = t106 * t103;
	t107 = t106 * t105;
	t102 = cos(pkin(6));
	t1 = [-t102 * t111 + t107, -t102 * t110 - t108, t112, t106 * pkin(1) + pkin(8) * t112 + 0; t102 * t108 + t110, t102 * t107 - t111, -t109, t104 * pkin(1) - pkin(8) * t109 + 0; t101 * t103, t101 * t105, t102, t102 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t116 = sin(pkin(6));
	t121 = cos(qJ(3));
	t131 = t116 * t121;
	t118 = sin(qJ(3));
	t130 = t118 * t116;
	t119 = sin(qJ(2));
	t129 = t119 * t121;
	t120 = sin(qJ(1));
	t128 = t120 * t119;
	t122 = cos(qJ(2));
	t127 = t120 * t122;
	t123 = cos(qJ(1));
	t126 = t123 * t119;
	t125 = t123 * t122;
	t124 = pkin(2) * t119 - pkin(9) * t122;
	t117 = cos(pkin(6));
	t115 = t122 * pkin(2) + t119 * pkin(9) + pkin(1);
	t114 = -t117 * t126 - t127;
	t113 = t116 * pkin(8) - t124 * t117;
	t1 = [(-t117 * t129 + t130) * t120 + t121 * t125, (t117 * t128 - t125) * t118 + t120 * t131, t117 * t127 + t126, t113 * t120 + t115 * t123 + 0; -t114 * t121 - t123 * t130, t114 * t118 - t123 * t131, -t117 * t125 + t128, -t113 * t123 + t115 * t120 + 0; t116 * t129 + t117 * t118, t117 * t121 - t119 * t130, -t116 * t122, t117 * pkin(8) + t124 * t116 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t140 = sin(pkin(6));
	t144 = sin(qJ(2));
	t156 = t140 * t144;
	t145 = sin(qJ(1));
	t155 = t140 * t145;
	t147 = cos(qJ(1));
	t154 = t140 * t147;
	t153 = t145 * t144;
	t146 = cos(qJ(2));
	t152 = t145 * t146;
	t151 = t147 * t144;
	t150 = t147 * t146;
	t149 = sin(qJ(3)) * pkin(3) + pkin(8);
	t136 = cos(qJ(3)) * pkin(3) + pkin(2);
	t142 = qJ(4) + pkin(9);
	t148 = t136 * t144 - t142 * t146;
	t141 = cos(pkin(6));
	t139 = qJ(3) + pkin(11);
	t138 = cos(t139);
	t137 = sin(t139);
	t135 = t141 * t151 + t152;
	t134 = t141 * t153 - t150;
	t133 = t136 * t146 + t142 * t144 + pkin(1);
	t132 = t140 * t149 - t148 * t141;
	t1 = [-t134 * t138 + t137 * t155, t134 * t137 + t138 * t155, t141 * t152 + t151, t132 * t145 + t133 * t147 + 0; t135 * t138 - t137 * t154, -t135 * t137 - t138 * t154, -t141 * t150 + t153, -t132 * t147 + t133 * t145 + 0; t141 * t137 + t138 * t156, -t137 * t156 + t141 * t138, -t140 * t146, t148 * t140 + t149 * t141 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (86->34), mult. (114->49), div. (0->0), fcn. (135->12), ass. (0->32)
	t170 = sin(pkin(11));
	t172 = cos(pkin(11));
	t165 = pkin(4) * t172 + qJ(5) * t170 + pkin(3);
	t175 = sin(qJ(3));
	t160 = t165 * t175 + pkin(8);
	t166 = -t170 * pkin(4) + qJ(5) * t172;
	t163 = t166 * t175 + pkin(2);
	t171 = sin(pkin(6));
	t173 = cos(pkin(6));
	t176 = sin(qJ(2));
	t178 = cos(qJ(3));
	t174 = qJ(4) + pkin(9);
	t179 = cos(qJ(2));
	t187 = t174 * t179;
	t191 = (t163 * t176 - t187) * t173 + (t173 * t165 * t176 + t171 * t166) * t178 - t171 * t160;
	t190 = t171 * t176;
	t177 = sin(qJ(1));
	t189 = t171 * t177;
	t180 = cos(qJ(1));
	t188 = t171 * t180;
	t186 = t177 * t176;
	t185 = t177 * t179;
	t184 = t180 * t176;
	t183 = t180 * t179;
	t169 = qJ(3) + pkin(11);
	t168 = cos(t169);
	t167 = sin(t169);
	t162 = t173 * t184 + t185;
	t161 = t173 * t186 - t183;
	t159 = t165 * t178 + t163;
	t157 = t159 * t179 + t174 * t176 + pkin(1);
	t1 = [t173 * t185 + t184, t161 * t168 - t167 * t189, -t161 * t167 - t168 * t189, t157 * t180 - t191 * t177 + 0; -t173 * t183 + t186, -t162 * t168 + t167 * t188, t162 * t167 + t168 * t188, t157 * t177 + t191 * t180 + 0; -t171 * t179, -t173 * t167 - t168 * t190, t167 * t190 - t173 * t168, pkin(7) + 0 + (t159 * t176 - t187) * t171 + (-t166 * t178 + t160) * t173; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:09
	% EndTime: 2020-11-04 22:24:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (125->40), mult. (152->61), div. (0->0), fcn. (186->14), ass. (0->41)
	t210 = sin(pkin(11));
	t212 = cos(pkin(11));
	t222 = pkin(4) + pkin(10);
	t203 = qJ(5) * t210 + t222 * t212 + pkin(3);
	t215 = sin(qJ(3));
	t196 = t203 * t215 + pkin(8);
	t204 = qJ(5) * t212 - t210 * t222;
	t197 = t204 * t215 + pkin(2);
	t211 = sin(pkin(6));
	t213 = cos(pkin(6));
	t216 = sin(qJ(2));
	t219 = cos(qJ(3));
	t208 = qJ(4) + pkin(5) + pkin(9);
	t220 = cos(qJ(2));
	t235 = t208 * t220;
	t236 = (t197 * t216 - t235) * t213 + (t203 * t213 * t216 + t211 * t204) * t219 - t211 * t196;
	t234 = t211 * t216;
	t217 = sin(qJ(1));
	t233 = t211 * t217;
	t232 = t211 * t220;
	t221 = cos(qJ(1));
	t231 = t211 * t221;
	t230 = t217 * t216;
	t229 = t217 * t220;
	t228 = t221 * t216;
	t227 = t221 * t220;
	t198 = t213 * t230 - t227;
	t209 = qJ(3) + pkin(11);
	t206 = sin(t209);
	t207 = cos(t209);
	t224 = -t198 * t206 - t207 * t233;
	t200 = t213 * t228 + t229;
	t223 = t200 * t206 + t207 * t231;
	t218 = cos(qJ(6));
	t214 = sin(qJ(6));
	t201 = t213 * t229 + t228;
	t199 = t213 * t227 - t230;
	t195 = t206 * t234 - t213 * t207;
	t194 = t203 * t219 + t197;
	t192 = t194 * t220 + t208 * t216 + pkin(1);
	t1 = [t201 * t218 + t224 * t214, -t201 * t214 + t224 * t218, -t198 * t207 + t206 * t233, t192 * t221 - t236 * t217 + 0; -t218 * t199 + t223 * t214, t214 * t199 + t223 * t218, t200 * t207 - t206 * t231, t192 * t217 + t236 * t221 + 0; t195 * t214 - t218 * t232, t195 * t218 + t214 * t232, t213 * t206 + t207 * t234, pkin(7) + 0 + (t194 * t216 - t235) * t211 + (-t204 * t219 + t196) * t213; 0, 0, 0, 1;];
	Tc_mdh = t1;
end