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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:19
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (107->0), mult. (416->0), div. (134->0), fcn. (472->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (460->23), mult. (1437->54), div. (151->8), fcn. (1823->9), ass. (0->33)
	t101 = sin(qJ(1));
	t79 = cos(qJ(1));
	t94 = cos(pkin(8));
	t76 = sin(pkin(8));
	t96 = t76 * cos(pkin(7));
	t72 = t101 * t96 + t79 * t94;
	t73 = -t101 * t94 + t79 * t96;
	t66 = atan2(t72, t73);
	t61 = sin(t66);
	t62 = cos(t66);
	t59 = t61 * t72 + t62 * t73;
	t58 = 0.1e1 / t59 ^ 2;
	t69 = t72 ^ 2;
	t70 = 0.1e1 / t73 ^ 2;
	t65 = t69 * t70 + 0.1e1;
	t63 = 0.1e1 / t65;
	t55 = t65 * t63;
	t104 = t55 - 0.1e1;
	t77 = sin(pkin(7));
	t103 = t77 ^ 2 * t76 ^ 2 * t58;
	t102 = 0.1e1 / t73;
	t68 = t73 * qJD(1);
	t98 = t70 * t72;
	t67 = t72 * qJD(1);
	t99 = t67 * t102 * t70;
	t100 = (t68 * t98 + t69 * t99) / t65 ^ 2;
	t54 = (t68 * t102 + t67 * t98) * t63;
	t91 = t54 * t104;
	t89 = -t61 * t73 + t62 * t72;
	t57 = 0.1e1 + t103;
	t56 = 0.1e1 / t57;
	t51 = -0.2e1 * t100 + 0.2e1 * (t63 * t68 * t70 + (-t70 * t100 + t63 * t99) * t72) * t72;
	t1 = [t51, 0, 0, 0, 0; (-((t51 * t72 + t55 * t68 - t73 * t91 - t68) * t62 + (-t51 * t73 + t55 * t67 - t72 * t91 - t67) * t61) * t56 + 0.2e1 * (t56 - 0.1e1 / t57 ^ 2 * t103) / t59 * t104 * t89 * (t89 * t54 + t61 * t68 - t62 * t67)) * t76 * t77 * t58, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (1004->38), mult. (3127->83), div. (198->12), fcn. (3985->13), ass. (0->54)
	t165 = sin(pkin(9));
	t166 = sin(pkin(8));
	t172 = sin(qJ(1));
	t174 = cos(qJ(1));
	t169 = cos(pkin(8));
	t170 = cos(pkin(7));
	t190 = t169 * t170;
	t167 = sin(pkin(7));
	t168 = cos(pkin(9));
	t192 = t167 * t168;
	t158 = -(t172 * t166 + t174 * t190) * t165 + t174 * t192;
	t151 = 0.1e1 / t158 ^ 2;
	t155 = t172 * t192 - (-t174 * t166 + t172 * t190) * t165;
	t204 = t151 * t155 ^ 2;
	t139 = 0.1e1 + t204;
	t137 = 0.1e1 / t139;
	t150 = 0.1e1 / t158;
	t186 = t150 * t158 + t204;
	t125 = t186 * t137;
	t206 = t125 - 0.1e1;
	t141 = t155 * qJD(1);
	t143 = t158 * qJD(1);
	t195 = t151 * t155;
	t124 = (t141 * t195 + t143 * t150) * t137;
	t205 = t206 * t124;
	t140 = atan2(-t155, -t158);
	t134 = sin(t140);
	t135 = cos(t140);
	t187 = t134 * t158 - t135 * t155;
	t131 = -t134 * t155 - t135 * t158;
	t130 = 0.1e1 / t131 ^ 2;
	t203 = t150 * t204;
	t191 = t167 * t169;
	t160 = t165 * t191 + t170 * t168;
	t202 = t160 ^ 2 * t130;
	t161 = -t170 * t165 + t168 * t191;
	t171 = sin(qJ(5));
	t173 = cos(qJ(5));
	t193 = t166 * t167;
	t154 = t161 * t173 + t171 * t193;
	t147 = 0.1e1 / t154 ^ 2;
	t200 = qJD(5) * t147;
	t185 = -t161 * t171 + t173 * t193;
	t146 = t185 ^ 2;
	t136 = t146 * t147 + 0.1e1;
	t196 = t154 * t200;
	t197 = t185 / t154 * t200;
	t199 = (-t146 * t197 - t185 * t196) / t136 ^ 2;
	t189 = t143 * t195;
	t132 = 0.1e1 / t136;
	t127 = 0.1e1 + t202;
	t126 = 0.1e1 / t127;
	t121 = -0.2e1 * t186 / t139 ^ 2 * (t141 * t203 + t189) + (0.2e1 * t189 - (-t151 * t158 + t150 - 0.2e1 * t203) * t141) * t137;
	t1 = [t121, 0, 0, 0, 0; (-((-t121 * t155 - t125 * t143 + t158 * t205 + t143) * t135 + (t121 * t158 - t125 * t141 + t155 * t205 + t141) * t134) * t126 + 0.2e1 * (t126 - 0.1e1 / t127 ^ 2 * t202) / t131 * t206 * t187 * (t187 * t124 - t134 * t143 + t135 * t141)) * t160 * t130, 0, 0, 0, 0; 0, 0, 0, 0, -0.2e1 * t199 - 0.2e1 * (t132 * t196 - (-t132 * t197 - t147 * t199) * t185) * t185;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end