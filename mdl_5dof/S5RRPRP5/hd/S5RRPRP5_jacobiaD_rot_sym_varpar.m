% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP5
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
%   Wie in S5RRPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:45:03
	% EndTime: 2019-12-29 18:45:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:45:08
	% EndTime: 2019-12-29 18:45:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:45:03
	% EndTime: 2019-12-29 18:45:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:45:03
	% EndTime: 2019-12-29 18:45:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:45:08
	% EndTime: 2019-12-29 18:45:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:44:59
	% EndTime: 2019-12-29 18:45:00
	% DurationCPUTime: 1.36s
	% Computational Cost: add. (4829->73), mult. (2860->163), div. (736->13), fcn. (3352->7), ass. (0->75)
	t121 = sin(qJ(1));
	t116 = t121 ^ 2;
	t113 = qJ(2) + pkin(8) + qJ(4);
	t111 = sin(t113);
	t107 = t111 ^ 2;
	t112 = cos(t113);
	t109 = 0.1e1 / t112 ^ 2;
	t155 = t107 * t109;
	t102 = t116 * t155 + 0.1e1;
	t106 = t111 * t107;
	t108 = 0.1e1 / t112;
	t110 = t108 * t109;
	t114 = qJD(2) + qJD(4);
	t154 = t108 * t111;
	t131 = t114 * (t106 * t110 + t154);
	t122 = cos(qJ(1));
	t146 = qJD(1) * t122;
	t138 = t121 * t146;
	t163 = 0.1e1 / t102 ^ 2 * (t116 * t131 + t138 * t155);
	t172 = -0.2e1 * t163;
	t100 = 0.1e1 / t102;
	t136 = 0.1e1 + t155;
	t170 = t121 * t136;
	t95 = t100 * t170;
	t171 = t121 * t95 - 0.1e1;
	t118 = 0.1e1 / t122;
	t149 = t118 * t121;
	t117 = t122 ^ 2;
	t169 = qJD(1) * (t116 / t117 + 0.1e1) * t149;
	t148 = t121 * t111;
	t99 = atan2(-t148, -t112);
	t97 = sin(t99);
	t142 = t97 * t148;
	t98 = cos(t99);
	t94 = -t112 * t98 - t142;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t152 = t112 * t114;
	t141 = t92 * t152;
	t151 = t114 * t121;
	t160 = t112 * t97;
	t139 = t109 * t151;
	t86 = (-(-t111 * t146 - t112 * t151) * t108 + t107 * t139) * t100;
	t81 = (t86 - t151) * t160 + (-t97 * t146 + (-t121 * t86 + t114) * t98) * t111;
	t167 = t81 * t91 * t92;
	t89 = t107 * t117 * t92 + 0.1e1;
	t168 = (t117 * t111 * t141 + (-t117 * t167 - t138 * t92) * t107) / t89 ^ 2;
	t87 = 0.1e1 / t89;
	t165 = t87 * t92;
	t164 = t91 * t87;
	t162 = t111 * t97;
	t161 = t111 * t98;
	t159 = t114 * t95;
	t157 = t122 * t92;
	t156 = t107 * t108;
	t153 = t111 * t122;
	t150 = t116 / t122 ^ 2;
	t147 = qJD(1) * t121;
	t145 = 0.2e1 * t167;
	t105 = t109 * t150 + 0.1e1;
	t144 = 0.2e1 / t105 ^ 2 * (t110 * t111 * t114 * t150 + t109 * t169);
	t143 = t91 * t168;
	t140 = t121 * t156;
	t137 = 0.2e1 * t92 * t168;
	t135 = 0.1e1 + t150;
	t134 = t136 * t122;
	t132 = t135 * t111 * t109;
	t130 = -t140 * t98 + t162;
	t103 = 0.1e1 / t105;
	t85 = (t100 * t130 - t162) * t122;
	t84 = t109 * t118 * t144 * t148 + ((-0.2e1 * t107 * t110 - t108) * t114 * t149 - qJD(1) * t132) * t103;
	t83 = (-t121 + t95) * t160 - t171 * t161;
	t82 = t170 * t172 + (qJD(1) * t134 + 0.2e1 * t121 * t131) * t100;
	t79 = (-t147 * t164 + (-0.2e1 * t143 + (-t114 * t83 - t81) * t165) * t122) * t112 + (t83 * t122 * t137 + (-t122 * t114 * t91 - ((-t121 * t82 - t146 * t95) * t98 + (t171 * t86 + t151 - t159) * t97) * t92 * t153 + (t122 * t145 + t147 * t92) * t83 - ((t82 - t146) * t97 + (t86 * t95 + t114 + (-t86 - t159) * t121) * t98) * t112 * t157) * t87) * t111;
	t1 = [t108 * t153 * t172 + (t114 * t134 - t147 * t154) * t100, t82, 0, t82, 0; (-t152 * t164 + (0.2e1 * t143 + (qJD(1) * t85 + t81) * t165) * t111) * t121 + (t85 * t137 * t111 + (-t85 * t141 + (t85 * t145 + (t86 * t161 + t97 * t152 + 0.2e1 * t130 * t163 + ((-t140 * t86 - t152) * t97 + (t106 * t139 - (t86 - 0.2e1 * t151) * t111) * t98) * t100) * t157) * t111 + (-t91 + (-t142 + (t142 - (t116 - t117) * t98 * t156) * t100) * t92) * t111 * qJD(1)) * t87) * t122, t79, 0, t79, 0; t135 * t108 * t144 + (-0.2e1 * t108 * t169 - t114 * t132) * t103, t84, 0, t84, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end