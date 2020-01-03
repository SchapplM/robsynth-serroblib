% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP2
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
%   Wie in S5RRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (118->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:07
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (3227->71), mult. (2645->159), div. (674->13), fcn. (3179->7), ass. (0->77)
	t106 = qJ(1) + qJ(2) + pkin(8);
	t104 = sin(t106);
	t140 = qJD(4) * t104;
	t105 = cos(t106);
	t114 = cos(qJ(4));
	t110 = 0.1e1 / t114;
	t139 = qJD(4) * t114;
	t113 = sin(qJ(4));
	t109 = t113 ^ 2;
	t111 = 0.1e1 / t114 ^ 2;
	t143 = t109 * t111;
	t107 = qJD(1) + qJD(2);
	t144 = t107 * t113;
	t99 = t104 ^ 2;
	t97 = t143 * t99 + 0.1e1;
	t95 = 0.1e1 / t97;
	t78 = (-(-t104 * t139 - t105 * t144) * t110 + t140 * t143) * t95;
	t164 = t78 - t140;
	t101 = 0.1e1 / t105;
	t148 = t101 * t104;
	t100 = t105 ^ 2;
	t163 = t107 * (0.1e1 / t100 * t99 + 0.1e1) * t148;
	t146 = t104 * t113;
	t94 = atan2(-t146, -t114);
	t92 = sin(t94);
	t134 = t92 * t146;
	t93 = cos(t94);
	t87 = -t114 * t93 - t134;
	t84 = 0.1e1 / t87;
	t85 = 0.1e1 / t87 ^ 2;
	t161 = 0.2e1 * t113;
	t160 = t95 - 0.1e1;
	t145 = t105 * t107;
	t126 = t104 * t109 * t145;
	t132 = t85 * t139;
	t151 = t114 * t92;
	t153 = t104 * t78;
	t73 = t164 * t151 + (-t92 * t145 + (qJD(4) - t153) * t93) * t113;
	t158 = t73 * t84 * t85;
	t82 = t100 * t109 * t85 + 0.1e1;
	t159 = (-t85 * t126 + (-t109 * t158 + t113 * t132) * t100) / t82 ^ 2;
	t80 = 0.1e1 / t82;
	t157 = t80 * t85;
	t108 = t113 * t109;
	t112 = t110 * t111;
	t142 = t110 * t113;
	t122 = qJD(4) * (t108 * t112 + t142);
	t156 = (t111 * t126 + t122 * t99) / t97 ^ 2;
	t155 = t84 * t80;
	t154 = 0.1e1 / t105 ^ 2 * t99;
	t152 = t105 * t85;
	t150 = t93 * t113;
	t129 = 0.1e1 + t143;
	t125 = t129 * t95;
	t88 = t104 * t125;
	t149 = qJD(4) * t88;
	t147 = t104 * t107;
	t141 = t111 * t113;
	t138 = 0.2e1 * t158;
	t91 = t111 * t154 + 0.1e1;
	t137 = 0.2e1 * (qJD(4) * t112 * t113 * t154 + t111 * t163) / t91 ^ 2;
	t136 = t84 * t159;
	t135 = t109 * t110 * t95;
	t133 = t113 * t160;
	t131 = 0.1e1 + t154;
	t130 = 0.2e1 * t85 * t159;
	t128 = -0.2e1 * t110 * t156;
	t127 = t104 * t135;
	t123 = t131 * t141;
	t89 = 0.1e1 / t91;
	t77 = (-t127 * t93 + t133 * t92) * t105;
	t76 = -t95 * t142 * t147 + (qJD(4) * t125 + t113 * t128) * t105;
	t75 = t88 * t151 + t150 + (-t150 * t88 - t151) * t104;
	t74 = t125 * t145 + 0.2e1 * (t122 * t95 - t129 * t156) * t104;
	t72 = -t89 * qJD(4) * t123 + (t131 * t137 - 0.2e1 * t163 * t89) * t110;
	t70 = (-t139 * t155 + (0.2e1 * t136 + (t107 * t77 + t73) * t157) * t113) * t104 + (t77 * t130 * t113 + (-t77 * t132 + (t77 * t138 + ((-t127 * t78 - t139 * t160 + t156 * t161) * t92 + (-t78 * t133 + (t109 * t128 + (t108 * t111 + t161) * t95 * qJD(4)) * t104) * t93) * t152) * t113 + (-t84 + (-(-t100 + t99) * t93 * t135 + t160 * t134) * t85) * t144) * t80) * t105;
	t1 = [t76, t76, 0, t74, 0; t70, t70, 0, (-t147 * t155 + (-0.2e1 * t136 + (-qJD(4) * t75 - t73) * t157) * t105) * t114 + (t75 * t105 * t130 + (-t105 * qJD(4) * t84 + (t105 * t138 + t147 * t85) * t75 + (-((t74 - t145) * t92 + (t78 * t88 + qJD(4) + (-t78 - t149) * t104) * t93) * t114 + (-(-t104 * t74 - t145 * t88) * t93 - (t153 * t88 - t149 - t164) * t92) * t113) * t152) * t80) * t113, 0; t72, t72, 0, t137 * t141 * t148 + (-t107 * t123 + (-0.2e1 * t109 * t112 - t110) * t101 * t140) * t89, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end