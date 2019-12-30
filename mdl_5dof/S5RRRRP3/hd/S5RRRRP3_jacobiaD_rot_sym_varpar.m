% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP3
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
%   Wie in S5RRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (195->0), mult. (111->0), div. (45->0), fcn. (66->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:02
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (4380->71), mult. (3455->159), div. (878->13), fcn. (4181->7), ass. (0->77)
	t109 = qJ(1) + qJ(2) + qJ(3);
	t106 = sin(t109);
	t142 = qJD(4) * t106;
	t107 = cos(t109);
	t116 = cos(qJ(4));
	t112 = 0.1e1 / t116;
	t141 = qJD(4) * t116;
	t115 = sin(qJ(4));
	t111 = t115 ^ 2;
	t113 = 0.1e1 / t116 ^ 2;
	t145 = t111 * t113;
	t108 = qJD(1) + qJD(2) + qJD(3);
	t146 = t108 * t115;
	t101 = t106 ^ 2;
	t99 = t101 * t145 + 0.1e1;
	t97 = 0.1e1 / t99;
	t80 = (-(-t106 * t141 - t107 * t146) * t112 + t142 * t145) * t97;
	t166 = t80 - t142;
	t103 = 0.1e1 / t107;
	t150 = t103 * t106;
	t102 = t107 ^ 2;
	t165 = t108 * (t101 / t102 + 0.1e1) * t150;
	t148 = t106 * t115;
	t96 = atan2(-t148, -t116);
	t93 = sin(t96);
	t136 = t93 * t148;
	t94 = cos(t96);
	t89 = -t116 * t94 - t136;
	t86 = 0.1e1 / t89;
	t87 = 0.1e1 / t89 ^ 2;
	t163 = 0.2e1 * t115;
	t162 = t97 - 0.1e1;
	t147 = t107 * t108;
	t128 = t106 * t111 * t147;
	t134 = t87 * t141;
	t154 = t116 * t93;
	t156 = t106 * t80;
	t75 = t166 * t154 + (-t93 * t147 + (qJD(4) - t156) * t94) * t115;
	t160 = t75 * t86 * t87;
	t84 = t102 * t111 * t87 + 0.1e1;
	t161 = (-t87 * t128 + (-t111 * t160 + t115 * t134) * t102) / t84 ^ 2;
	t82 = 0.1e1 / t84;
	t159 = t82 * t87;
	t110 = t115 * t111;
	t114 = t112 * t113;
	t144 = t112 * t115;
	t124 = qJD(4) * (t110 * t114 + t144);
	t158 = (t101 * t124 + t113 * t128) / t99 ^ 2;
	t157 = t86 * t82;
	t155 = t107 * t87;
	t153 = t94 * t115;
	t131 = 0.1e1 + t145;
	t127 = t131 * t97;
	t90 = t106 * t127;
	t152 = qJD(4) * t90;
	t151 = t101 / t107 ^ 2;
	t149 = t106 * t108;
	t143 = t113 * t115;
	t140 = 0.2e1 * t160;
	t95 = t113 * t151 + 0.1e1;
	t139 = 0.2e1 * (qJD(4) * t114 * t115 * t151 + t113 * t165) / t95 ^ 2;
	t138 = t86 * t161;
	t137 = t111 * t112 * t97;
	t135 = t162 * t115;
	t133 = 0.2e1 * t87 * t161;
	t132 = 0.1e1 + t151;
	t130 = -0.2e1 * t112 * t158;
	t129 = t106 * t137;
	t125 = t132 * t143;
	t91 = 0.1e1 / t95;
	t79 = (-t129 * t94 + t135 * t93) * t107;
	t78 = -t97 * t144 * t149 + (qJD(4) * t127 + t115 * t130) * t107;
	t77 = t90 * t154 + t153 + (-t153 * t90 - t154) * t106;
	t76 = t127 * t147 + 0.2e1 * (t124 * t97 - t131 * t158) * t106;
	t74 = -t91 * qJD(4) * t125 + (t132 * t139 - 0.2e1 * t165 * t91) * t112;
	t72 = (-t141 * t157 + (0.2e1 * t138 + (t108 * t79 + t75) * t159) * t115) * t106 + (t79 * t133 * t115 + (-t79 * t134 + (t79 * t140 + ((-t80 * t129 - t141 * t162 + t158 * t163) * t93 + (-t80 * t135 + (t111 * t130 + (t110 * t113 + t163) * t97 * qJD(4)) * t106) * t94) * t155) * t115 + (-t86 + (-(t101 - t102) * t94 * t137 + t162 * t136) * t87) * t146) * t82) * t107;
	t1 = [t78, t78, t78, t76, 0; t72, t72, t72, (-t149 * t157 + (-0.2e1 * t138 + (-qJD(4) * t77 - t75) * t159) * t107) * t116 + (t77 * t107 * t133 + (-t107 * qJD(4) * t86 + (t107 * t140 + t149 * t87) * t77 + (-((t76 - t147) * t93 + (t80 * t90 + qJD(4) + (-t80 - t152) * t106) * t94) * t116 + (-(-t106 * t76 - t147 * t90) * t94 - (t156 * t90 - t152 - t166) * t93) * t115) * t155) * t82) * t115, 0; t74, t74, t74, t139 * t143 * t150 + (-t108 * t125 + (-0.2e1 * t111 * t114 - t112) * t103 * t142) * t91, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end