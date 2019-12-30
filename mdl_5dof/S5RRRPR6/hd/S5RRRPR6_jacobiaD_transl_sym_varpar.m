% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR6_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:20
	% EndTime: 2019-12-29 20:02:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:20
	% EndTime: 2019-12-29 20:02:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:26
	% EndTime: 2019-12-29 20:02:26
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:18
	% EndTime: 2019-12-29 20:02:18
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(7) + pkin(6);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:22
	% EndTime: 2019-12-29 20:02:23
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (179->35), mult. (217->48), div. (0->0), fcn. (148->6), ass. (0->34)
	t183 = qJ(2) + qJ(3);
	t181 = cos(t183);
	t211 = r_i_i_C(3) + qJ(4);
	t195 = t211 * t181;
	t180 = sin(t183);
	t178 = t180 * qJD(4);
	t182 = qJD(2) + qJD(3);
	t214 = pkin(3) + r_i_i_C(1);
	t203 = t214 * t180;
	t184 = sin(qJ(2));
	t210 = pkin(2) * qJD(2);
	t204 = t184 * t210;
	t217 = (-t203 + t195) * t182 + (r_i_i_C(2) + pkin(7) + pkin(6)) * qJD(1) + t178 - t204;
	t213 = pkin(2) * t184;
	t209 = t181 * t182;
	t187 = cos(qJ(1));
	t208 = t182 * t187;
	t185 = sin(qJ(1));
	t207 = qJD(1) * t185;
	t206 = qJD(1) * t187;
	t205 = qJD(4) * t181;
	t202 = t214 * t187;
	t201 = t185 * t209;
	t200 = t185 * t205 + t206 * t195;
	t197 = t180 * t207;
	t199 = t187 * t205 + t214 * t197;
	t196 = t211 * t180;
	t194 = t211 * t185;
	t192 = -t214 * t181 - t196;
	t191 = -t182 * t203 + t211 * t209 + t178;
	t186 = cos(qJ(2));
	t190 = qJD(1) * (-t186 * pkin(2) - pkin(1) + t192);
	t189 = t192 * t182 - t186 * t210;
	t1 = [-t217 * t185 + t187 * t190, (-t195 + t213) * t207 + t189 * t187 + t199, -t196 * t208 + (-qJD(1) * t194 - t182 * t202) * t181 + t199, t181 * t208 - t197, 0; t185 * t190 + t217 * t187, (-t203 - t213) * t206 + t189 * t185 + t200, -t214 * t201 + (-qJD(1) * t202 - t182 * t194) * t180 + t200, t180 * t206 + t201, 0; 0, t191 - t204, t191, t182 * t180, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:21
	% EndTime: 2019-12-29 20:02:21
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (426->48), mult. (539->67), div. (0->0), fcn. (460->8), ass. (0->46)
	t100 = qJ(2) + qJ(3);
	t98 = cos(t100);
	t135 = qJ(4) * t98;
	t139 = pkin(3) + pkin(4);
	t97 = sin(t100);
	t99 = qJD(2) + qJD(3);
	t150 = t99 * t135 + (-t139 * t99 + qJD(4)) * t97;
	t142 = qJD(5) - t99;
	t101 = sin(qJ(5));
	t104 = cos(qJ(5));
	t148 = -t101 * t98 + t104 * t97;
	t140 = t142 * t148;
	t102 = sin(qJ(2));
	t136 = pkin(2) * qJD(2);
	t125 = t102 * t136;
	t149 = -(pkin(8) + r_i_i_C(3) - pkin(7) - pkin(6)) * qJD(1) - t125 + t150;
	t116 = t101 * t97 + t104 * t98;
	t84 = t142 * t116;
	t103 = sin(qJ(1));
	t106 = cos(qJ(1));
	t114 = qJD(1) * t148;
	t77 = -t103 * t114 - t84 * t106;
	t113 = qJD(1) * t116;
	t78 = t103 * t113 - t140 * t106;
	t145 = t77 * r_i_i_C(1) + t78 * r_i_i_C(2);
	t79 = t103 * t84 - t106 * t114;
	t80 = t140 * t103 + t106 * t113;
	t144 = -t79 * r_i_i_C(1) - t80 * r_i_i_C(2);
	t143 = -r_i_i_C(1) * t140 + t84 * r_i_i_C(2);
	t138 = pkin(2) * t102;
	t137 = t98 * t99;
	t132 = qJD(4) * t98;
	t131 = qJD(1) * t103;
	t130 = qJD(1) * t106;
	t124 = qJD(1) * t135;
	t123 = t97 * t131;
	t122 = t97 * t130;
	t119 = t103 * t132 + t106 * t124 - t144;
	t118 = t106 * t132 + t139 * t123 - t145;
	t115 = -qJ(4) * t97 - t139 * t98;
	t111 = t115 * t99;
	t110 = -t143 + t150;
	t105 = cos(qJ(2));
	t109 = qJD(1) * (-t105 * pkin(2) - pkin(1) + t115);
	t108 = -t105 * t136 + t111;
	t1 = [-t80 * r_i_i_C(1) + t79 * r_i_i_C(2) - t149 * t103 + t106 * t109, (-t135 + t138) * t131 + t108 * t106 + t118, -t103 * t124 + t106 * t111 + t118, t106 * t137 - t123, t145; -t78 * r_i_i_C(1) + t77 * r_i_i_C(2) + t103 * t109 + t149 * t106, (-t139 * t97 - t138) * t130 + t108 * t103 + t119, t103 * t111 - t139 * t122 + t119, t103 * t137 + t122, t144; 0, t110 - t125, t110, t99 * t97, t143;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end