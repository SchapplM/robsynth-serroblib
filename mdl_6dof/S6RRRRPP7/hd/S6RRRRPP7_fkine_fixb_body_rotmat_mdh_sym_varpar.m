% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t104 = cos(qJ(1));
	t103 = sin(qJ(1));
	t1 = [t104, -t103, 0, 0; t103, t104, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t105 = sin(pkin(6));
	t108 = sin(qJ(1));
	t116 = t108 * t105;
	t107 = sin(qJ(2));
	t115 = t108 * t107;
	t109 = cos(qJ(2));
	t114 = t108 * t109;
	t110 = cos(qJ(1));
	t113 = t110 * t105;
	t112 = t110 * t107;
	t111 = t110 * t109;
	t106 = cos(pkin(6));
	t1 = [-t106 * t115 + t111, -t106 * t114 - t112, t116, t110 * pkin(1) + pkin(8) * t116 + 0; t106 * t112 + t114, t106 * t111 - t115, -t113, t108 * pkin(1) - pkin(8) * t113 + 0; t105 * t107, t105 * t109, t106, t106 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t120 = sin(pkin(6));
	t125 = cos(qJ(3));
	t135 = t120 * t125;
	t122 = sin(qJ(3));
	t134 = t122 * t120;
	t123 = sin(qJ(2));
	t133 = t123 * t125;
	t124 = sin(qJ(1));
	t132 = t124 * t123;
	t126 = cos(qJ(2));
	t131 = t124 * t126;
	t127 = cos(qJ(1));
	t130 = t127 * t123;
	t129 = t127 * t126;
	t128 = pkin(2) * t123 - pkin(9) * t126;
	t121 = cos(pkin(6));
	t119 = t126 * pkin(2) + t123 * pkin(9) + pkin(1);
	t118 = t121 * t130 + t131;
	t117 = t120 * pkin(8) - t128 * t121;
	t1 = [(-t121 * t133 + t134) * t124 + t125 * t129, (t121 * t132 - t129) * t122 + t124 * t135, t121 * t131 + t130, t117 * t124 + t119 * t127 + 0; t118 * t125 - t127 * t134, -t118 * t122 - t127 * t135, -t121 * t129 + t132, -t117 * t127 + t119 * t124 + 0; t120 * t133 + t121 * t122, t121 * t125 - t123 * t134, -t120 * t126, t121 * pkin(8) + t128 * t120 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t143 = sin(pkin(6));
	t150 = cos(qJ(3));
	t164 = t143 * t150;
	t145 = sin(qJ(4));
	t151 = cos(qJ(2));
	t163 = t145 * t151;
	t146 = sin(qJ(3));
	t162 = t146 * t143;
	t147 = sin(qJ(2));
	t161 = t147 * t150;
	t148 = sin(qJ(1));
	t160 = t148 * t151;
	t149 = cos(qJ(4));
	t159 = t149 * t151;
	t158 = t150 * t151;
	t152 = cos(qJ(1));
	t157 = t152 * t147;
	t156 = t152 * t151;
	t141 = t150 * pkin(3) + t146 * pkin(10) + pkin(2);
	t155 = t151 * pkin(9) - t141 * t147;
	t154 = t146 * pkin(3) - t150 * pkin(10) + pkin(8);
	t144 = cos(pkin(6));
	t139 = t144 * t161 - t162;
	t153 = t139 * t145 + t144 * t159;
	t140 = t145 * t158 - t149 * t147;
	t138 = t143 * t161 + t144 * t146;
	t137 = t147 * pkin(9) + t141 * t151 + pkin(1);
	t136 = t143 * t154 + t155 * t144;
	t1 = [(-t139 * t148 + t150 * t156) * t149 + (t144 * t160 + t157) * t145, -t152 * t140 + t153 * t148, (-t148 * t144 * t147 + t156) * t146 - t148 * t164, t136 * t148 + t137 * t152 + 0; (t139 * t149 - t144 * t163) * t152 + t148 * (t147 * t145 + t149 * t158), -t148 * t140 - t153 * t152, (t144 * t157 + t160) * t146 + t152 * t164, -t136 * t152 + t137 * t148 + 0; t138 * t149 - t143 * t163, -t138 * t145 - t143 * t159, -t144 * t150 + t147 * t162, -t155 * t143 + t154 * t144 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:46
	% EndTime: 2020-11-04 22:36:47
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t178 = sin(pkin(6));
	t184 = cos(qJ(3));
	t198 = t178 * t184;
	t185 = cos(qJ(2));
	t197 = t178 * t185;
	t181 = sin(qJ(3));
	t196 = t181 * t178;
	t182 = sin(qJ(2));
	t195 = t182 * t184;
	t183 = sin(qJ(1));
	t194 = t183 * t182;
	t193 = t183 * t185;
	t186 = cos(qJ(1));
	t192 = t186 * t182;
	t191 = t186 * t185;
	t174 = cos(qJ(4)) * pkin(4) + pkin(3);
	t180 = qJ(5) + pkin(10);
	t167 = t174 * t184 + t180 * t181 + pkin(2);
	t173 = sin(qJ(4)) * pkin(4) + pkin(9);
	t190 = t167 * t182 - t173 * t185;
	t189 = t174 * t181 - t180 * t184 + pkin(8);
	t179 = cos(pkin(6));
	t169 = t179 * t195 - t196;
	t188 = -t169 * t183 + t184 * t191;
	t187 = t169 * t186 + t184 * t193;
	t177 = qJ(4) + pkin(11);
	t176 = cos(t177);
	t175 = sin(t177);
	t171 = t179 * t193 + t192;
	t170 = t179 * t191 - t194;
	t168 = t178 * t195 + t179 * t181;
	t166 = t167 * t185 + t173 * t182 + pkin(1);
	t165 = t178 * t189 - t190 * t179;
	t1 = [t171 * t175 + t188 * t176, t171 * t176 - t188 * t175, (-t179 * t194 + t191) * t181 - t183 * t198, t165 * t183 + t166 * t186 + 0; -t175 * t170 + t187 * t176, -t176 * t170 - t187 * t175, (t179 * t192 + t193) * t181 + t186 * t198, -t165 * t186 + t166 * t183 + 0; t168 * t176 - t175 * t197, -t168 * t175 - t176 * t197, -t179 * t184 + t182 * t196, t190 * t178 + t189 * t179 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:47
	% EndTime: 2020-11-04 22:36:47
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (142->41), mult. (199->61), div. (0->0), fcn. (233->14), ass. (0->40)
	t218 = qJ(5) + pkin(10);
	t220 = sin(qJ(3));
	t224 = cos(qJ(3));
	t214 = sin(pkin(11));
	t216 = cos(pkin(11));
	t209 = pkin(5) * t216 + qJ(6) * t214 + pkin(4);
	t210 = -t214 * pkin(5) + qJ(6) * t216;
	t219 = sin(qJ(4));
	t223 = cos(qJ(4));
	t240 = -t209 * t223 - t210 * t219 - pkin(3);
	t241 = t218 * t224 + t240 * t220 - pkin(8);
	t215 = sin(pkin(6));
	t239 = t215 * t224;
	t225 = cos(qJ(2));
	t238 = t215 * t225;
	t237 = t220 * t215;
	t221 = sin(qJ(2));
	t236 = t221 * t224;
	t222 = sin(qJ(1));
	t235 = t222 * t221;
	t234 = t222 * t225;
	t226 = cos(qJ(1));
	t233 = t226 * t221;
	t232 = t226 * t225;
	t228 = t209 * t219 - t210 * t223 + pkin(9);
	t201 = t218 * t220 - t224 * t240 + pkin(2);
	t227 = t201 * t221 - t228 * t225;
	t217 = cos(pkin(6));
	t213 = qJ(4) + pkin(11);
	t212 = cos(t213);
	t211 = sin(t213);
	t207 = t217 * t234 + t233;
	t206 = t217 * t232 - t235;
	t205 = t217 * t236 - t237;
	t204 = t215 * t236 + t217 * t220;
	t203 = -t205 * t222 + t224 * t232;
	t202 = t205 * t226 + t224 * t234;
	t200 = t201 * t225 + t228 * t221 + pkin(1);
	t199 = t215 * t241 + t227 * t217;
	t1 = [t203 * t212 + t207 * t211, (-t217 * t235 + t232) * t220 - t222 * t239, t203 * t211 - t207 * t212, -t199 * t222 + t200 * t226 + 0; t202 * t212 - t211 * t206, (t217 * t233 + t234) * t220 + t226 * t239, t202 * t211 + t212 * t206, t199 * t226 + t200 * t222 + 0; t204 * t212 - t211 * t238, -t217 * t224 + t221 * t237, t204 * t211 + t212 * t238, t227 * t215 - t241 * t217 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end