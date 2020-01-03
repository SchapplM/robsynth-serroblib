% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.06s
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
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-12-31 20:59:09
	% EndTime: 2019-12-31 20:59:09
	% DurationCPUTime: 0.16s
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
	t198 = t180 * t207;
	t199 = t187 * t205 + t214 * t198;
	t196 = t211 * t180;
	t194 = t211 * t185;
	t192 = -t214 * t181 - t196;
	t191 = -t182 * t203 + t211 * t209 + t178;
	t186 = cos(qJ(2));
	t190 = qJD(1) * (-pkin(2) * t186 - pkin(1) + t192);
	t189 = t192 * t182 - t186 * t210;
	t1 = [-t217 * t185 + t187 * t190, (-t195 + t213) * t207 + t189 * t187 + t199, -t196 * t208 + (-qJD(1) * t194 - t182 * t202) * t181 + t199, t181 * t208 - t198, 0; t185 * t190 + t217 * t187, (-t203 - t213) * t206 + t189 * t185 + t200, -t214 * t201 + (-qJD(1) * t202 - t182 * t194) * t180 + t200, t180 * t206 + t201, 0; 0, t191 - t204, t191, t182 * t180, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (221->39), mult. (265->49), div. (0->0), fcn. (180->6), ass. (0->34)
	t62 = qJ(2) + qJ(3);
	t59 = sin(t62);
	t57 = t59 * qJD(4);
	t61 = qJD(2) + qJD(3);
	t85 = pkin(3) + pkin(4) + r_i_i_C(1);
	t75 = t85 * t59;
	t63 = sin(qJ(2));
	t89 = pkin(2) * qJD(2);
	t82 = t63 * t89;
	t60 = cos(t62);
	t90 = r_i_i_C(2) + qJ(4);
	t98 = t90 * t60;
	t99 = (-t98 + t75) * t61 + (r_i_i_C(3) + qJ(5) - pkin(7) - pkin(6)) * qJD(1) - t57 + t82;
	t94 = pkin(2) * t63;
	t64 = sin(qJ(1));
	t92 = t60 * t64;
	t66 = cos(qJ(1));
	t91 = t61 * t66;
	t88 = qJD(1) * t64;
	t87 = qJD(1) * t66;
	t86 = qJD(4) * t60;
	t83 = t64 * t86 + t87 * t98;
	t80 = t59 * t88;
	t79 = t90 * t59;
	t77 = t90 * t64;
	t76 = t66 * t86 + t85 * t80;
	t74 = t85 * t61;
	t73 = t85 * t66;
	t71 = -t59 * t74 + t61 * t98 + t57;
	t70 = -t60 * t85 - t79;
	t65 = cos(qJ(2));
	t69 = -qJD(5) + (-t65 * pkin(2) - pkin(1) + t70) * qJD(1);
	t68 = t61 * t70 - t65 * t89;
	t1 = [t99 * t64 + t69 * t66, (-t98 + t94) * t88 + t68 * t66 + t76, -t79 * t91 + (-qJD(1) * t77 - t61 * t73) * t60 + t76, t60 * t91 - t80, -t87; t69 * t64 - t99 * t66, (-t75 - t94) * t87 + t68 * t64 + t83, -t74 * t92 + (-qJD(1) * t73 - t61 * t77) * t59 + t83, t59 * t87 + t61 * t92, -t88; 0, t71 - t82, t71, t61 * t59, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end