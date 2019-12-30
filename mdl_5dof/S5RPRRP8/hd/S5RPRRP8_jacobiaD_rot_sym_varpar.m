% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RPRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:29
	% EndTime: 2019-12-29 17:24:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:29
	% EndTime: 2019-12-29 17:24:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:29
	% EndTime: 2019-12-29 17:24:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:24
	% EndTime: 2019-12-29 17:24:24
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (128->13), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->21)
	t54 = sin(qJ(3));
	t55 = cos(qJ(1));
	t67 = sin(qJ(1));
	t68 = cos(qJ(3));
	t59 = t67 * t54 + t55 * t68;
	t44 = 0.1e1 / t59 ^ 2;
	t70 = t44 * t59;
	t69 = qJD(1) - qJD(3);
	t43 = 0.1e1 / t59;
	t58 = -t55 * t54 + t67 * t68;
	t38 = t69 * t58;
	t42 = t58 ^ 2;
	t64 = t42 * t44;
	t41 = 0.1e1 + t64;
	t65 = t69 * t70;
	t62 = t58 * t65;
	t45 = t43 * t44;
	t63 = t42 * t45;
	t66 = (t38 * t63 + t62) / t41 ^ 2;
	t39 = 0.1e1 / t41;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0.2e1 * (t43 * t59 + t64) * t66 + (-0.2e1 * t62 - (-t43 + 0.2e1 * t63 + t70) * t38) * t39, 0, -0.2e1 * t66 - 0.2e1 * (-t39 * t65 - (t38 * t39 * t45 - t44 * t66) * t58) * t58, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:24
	% EndTime: 2019-12-29 17:24:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:24:30
	% EndTime: 2019-12-29 17:24:32
	% DurationCPUTime: 2.80s
	% Computational Cost: add. (2108->83), mult. (5286->193), div. (692->13), fcn. (6766->9), ass. (0->85)
	t213 = sin(qJ(3));
	t214 = sin(qJ(1));
	t215 = cos(qJ(3));
	t216 = cos(qJ(1));
	t129 = t216 * t213 - t214 * t215;
	t123 = t129 ^ 2;
	t144 = sin(qJ(4));
	t140 = t144 ^ 2;
	t145 = cos(qJ(4));
	t142 = 0.1e1 / t145 ^ 2;
	t192 = t140 * t142;
	t122 = t123 * t192 + 0.1e1;
	t120 = 0.1e1 / t122;
	t141 = 0.1e1 / t145;
	t185 = qJD(4) * t145;
	t173 = t129 * t185;
	t187 = qJD(4) * t129;
	t176 = t142 * t187;
	t128 = -t214 * t213 - t216 * t215;
	t228 = qJD(1) - qJD(3);
	t114 = t228 * t128;
	t204 = t114 * t144;
	t102 = ((t173 + t204) * t141 + t140 * t176) * t120;
	t194 = t129 * t144;
	t119 = atan2(t194, t145);
	t117 = sin(t119);
	t118 = cos(t119);
	t111 = t117 * t194 + t118 * t145;
	t109 = 0.1e1 / t111 ^ 2;
	t124 = t128 ^ 2;
	t196 = t124 * t140;
	t106 = t109 * t196 + 0.1e1;
	t104 = 0.1e1 / t106;
	t108 = 0.1e1 / t111;
	t115 = t228 * t129;
	t139 = t144 * t140;
	t193 = t140 * t141;
	t178 = t129 * t193;
	t167 = t118 * t178;
	t143 = t141 * t142;
	t191 = t141 * t144;
	t155 = qJD(4) * (t139 * t143 + t191);
	t180 = t114 * t129 * t142;
	t207 = (t123 * t155 + t140 * t180) / t122 ^ 2;
	t183 = 0.2e1 * t207;
	t198 = t118 * t144;
	t97 = (t102 * t129 - qJD(4)) * t198 + (t204 + (-t102 + t187) * t145) * t117;
	t210 = t108 * t109 * t97;
	t184 = -0.2e1 * t210;
	t200 = t117 * t144;
	t230 = (-t200 + (-t167 + t200) * t120) * t128;
	t240 = t109 * t230;
	t212 = 0.1e1 / t106 ^ 2 * (-t196 * t210 + (-t115 * t128 * t140 + t124 * t144 * t185) * t109);
	t246 = -0.2e1 * t144 * t212;
	t247 = ((((((-t115 * t178 + (t114 * t193 + t139 * t176 + (-t102 + 0.2e1 * t187) * t144) * t128) * t118 - ((t102 * t178 + t185) * t128 - t115 * t144) * t117) * t120 + ((t102 * t118 + t117 * t183) * t128 - t115 * t117) * t144 - (-t117 * t185 + t167 * t183) * t128) * t109 - t230 * t184) * t144 - t185 * t240) * t104 - t240 * t246) * t128 + t129 * t108 * t246 + (t108 * t173 + (t108 * t114 + (t115 * t230 - t129 * t97) * t109) * t144) * t104;
	t237 = 0.2e1 * t141;
	t236 = -0.2e1 * t207;
	t171 = 0.1e1 + t192;
	t227 = t128 * t191 * t236 + (t171 * t128 * qJD(4) - t115 * t191) * t120;
	t224 = t129 * t171;
	t125 = 0.1e1 / t128;
	t107 = t120 * t224;
	t199 = t117 * t145;
	t98 = t129 * t199 - t198 + (t118 * t194 - t199) * t107;
	t211 = t109 * t98;
	t126 = 0.1e1 / t128 ^ 2;
	t197 = t123 * t126;
	t116 = t142 * t197 + 0.1e1;
	t186 = qJD(4) * t144;
	t127 = t125 / t124;
	t202 = t115 * t127;
	t208 = (t126 * t180 + (t126 * t143 * t186 + t142 * t202) * t123) / t116 ^ 2;
	t205 = t109 * t128;
	t195 = t126 * t129;
	t190 = t142 * t144;
	t189 = t107 - t129;
	t181 = t141 * t208;
	t177 = t129 * t190;
	t175 = t142 * t186;
	t172 = t107 * t129 - 0.1e1;
	t170 = -0.2e1 * t181;
	t157 = t125 * t128 + t197;
	t112 = 0.1e1 / t116;
	t96 = t224 * t236 + (t171 * t114 + 0.2e1 * t129 * t155) * t120;
	t1 = [t227, 0, -t227, t96, 0; t247, 0, -t247, 0.2e1 * (t108 * t145 - t144 * t211) * t128 * t212 + ((t115 * t108 + (qJD(4) * t98 + t97) * t205) * t145 + (-t115 * t211 + (qJD(4) * t108 + t98 * t184 + ((t107 * t114 + t129 * t96) * t198 + (t189 * qJD(4) - t172 * t102) * t200) * t109) * t128 + ((t114 - t96) * t117 + (t172 * qJD(4) - t189 * t102) * t118) * t145 * t205) * t144) * t104, 0; t112 * t175 + t170 + (t170 * t195 + (t114 * t126 * t237 + (t126 * t175 + t202 * t237) * t129) * t112) * t129, 0, 0.2e1 * t157 * t181 + (-t157 * t175 + (-0.2e1 * t114 * t195 + (-0.2e1 * t123 * t127 - t126 * t128 + t125) * t115) * t141) * t112, -0.2e1 * t125 * t177 * t208 + (t115 * t126 * t177 - (-t114 * t190 + (-0.2e1 * t140 * t143 - t141) * t187) * t125) * t112, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end