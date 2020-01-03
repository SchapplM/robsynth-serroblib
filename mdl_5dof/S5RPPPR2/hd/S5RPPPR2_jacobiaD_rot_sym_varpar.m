% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   Wie in S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (183->0), mult. (416->0), div. (134->0), fcn. (472->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (460->25), mult. (1437->54), div. (151->8), fcn. (1823->9), ass. (0->33)
	t101 = cos(pkin(8));
	t83 = sin(pkin(8));
	t103 = t83 * cos(pkin(7));
	t108 = sin(qJ(1));
	t86 = cos(qJ(1));
	t77 = t86 * t101 + t108 * t103;
	t79 = -t108 * t101 + t86 * t103;
	t71 = atan2(-t77, -t79);
	t66 = sin(t71);
	t67 = cos(t71);
	t64 = -t66 * t77 - t67 * t79;
	t63 = 0.1e1 / t64 ^ 2;
	t74 = t77 ^ 2;
	t75 = 0.1e1 / t79 ^ 2;
	t70 = t74 * t75 + 0.1e1;
	t68 = 0.1e1 / t70;
	t60 = t70 * t68;
	t111 = t60 - 0.1e1;
	t84 = sin(pkin(7));
	t110 = t84 ^ 2 * t83 ^ 2 * t63;
	t109 = 0.1e1 / t79;
	t105 = t75 * t77;
	t72 = t77 * qJD(1);
	t106 = t72 * t109 * t75;
	t73 = t79 * qJD(1);
	t107 = (t73 * t105 + t74 * t106) / t70 ^ 2;
	t59 = (t72 * t105 + t73 * t109) * t68;
	t98 = t59 * t111;
	t96 = t66 * t79 - t67 * t77;
	t62 = 0.1e1 + t110;
	t61 = 0.1e1 / t62;
	t56 = -0.2e1 * t107 + 0.2e1 * (t68 * t73 * t75 + (t68 * t106 - t75 * t107) * t77) * t77;
	t1 = [t56, 0, 0, 0, 0; (-((-t56 * t77 - t60 * t73 + t79 * t98 + t73) * t67 + (t56 * t79 - t60 * t72 + t77 * t98 + t72) * t66) * t61 + 0.2e1 * (t61 - 0.1e1 / t62 ^ 2 * t110) / t64 * t111 * t96 * (t96 * t59 - t66 * t73 + t67 * t72)) * t63 * t83 * t84, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (1004->38), mult. (3127->83), div. (198->12), fcn. (3985->13), ass. (0->53)
	t168 = sin(pkin(9));
	t169 = sin(pkin(8));
	t174 = cos(qJ(1));
	t200 = sin(pkin(7));
	t201 = cos(pkin(9));
	t189 = t200 * t201;
	t170 = cos(pkin(8));
	t171 = cos(pkin(7));
	t192 = t170 * t171;
	t202 = sin(qJ(1));
	t158 = -t174 * t189 + (t202 * t169 + t174 * t192) * t168;
	t152 = 0.1e1 / t158 ^ 2;
	t156 = -t202 * t189 + (-t174 * t169 + t202 * t192) * t168;
	t206 = t152 * t156 ^ 2;
	t140 = 0.1e1 + t206;
	t138 = 0.1e1 / t140;
	t151 = 0.1e1 / t158;
	t187 = t151 * t158 + t206;
	t126 = t187 * t138;
	t209 = t126 - 0.1e1;
	t142 = t156 * qJD(1);
	t144 = t158 * qJD(1);
	t195 = t152 * t156;
	t125 = (t142 * t195 + t144 * t151) * t138;
	t208 = t209 * t125;
	t141 = atan2(-t156, -t158);
	t135 = sin(t141);
	t136 = cos(t141);
	t188 = t135 * t158 - t136 * t156;
	t132 = -t135 * t156 - t136 * t158;
	t131 = 0.1e1 / t132 ^ 2;
	t161 = t200 * t170 * t168 + t171 * t201;
	t207 = t131 * t161 ^ 2;
	t205 = t151 * t206;
	t162 = -t171 * t168 + t170 * t189;
	t172 = sin(qJ(5));
	t173 = cos(qJ(5));
	t190 = t169 * t200;
	t155 = t162 * t173 + t172 * t190;
	t148 = 0.1e1 / t155 ^ 2;
	t203 = qJD(5) * t148;
	t184 = -t162 * t172 + t173 * t190;
	t147 = t184 ^ 2;
	t137 = t147 * t148 + 0.1e1;
	t196 = t155 * t203;
	t197 = t184 / t155 * t203;
	t199 = 0.1e1 / t137 ^ 2 * (-t147 * t197 - t184 * t196);
	t191 = t144 * t195;
	t133 = 0.1e1 / t137;
	t128 = 0.1e1 + t207;
	t127 = 0.1e1 / t128;
	t122 = -0.2e1 * t187 / t140 ^ 2 * (t142 * t205 + t191) + (0.2e1 * t191 + (t152 * t158 - t151 + 0.2e1 * t205) * t142) * t138;
	t1 = [t122, 0, 0, 0, 0; (-((-t122 * t156 - t126 * t144 + t158 * t208 + t144) * t136 + (t122 * t158 - t126 * t142 + t156 * t208 + t142) * t135) * t127 + 0.2e1 * (t127 - 0.1e1 / t128 ^ 2 * t207) / t132 * t209 * t188 * (t188 * t125 - t135 * t144 + t136 * t142)) * t161 * t131, 0, 0, 0, 0; 0, 0, 0, 0, -0.2e1 * t199 - 0.2e1 * (t133 * t196 - (-t133 * t197 - t148 * t199) * t184) * t184;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end