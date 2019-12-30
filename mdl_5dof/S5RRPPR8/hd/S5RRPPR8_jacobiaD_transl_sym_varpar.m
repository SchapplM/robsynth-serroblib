% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:25:54
	% EndTime: 2019-12-29 18:25:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:25:49
	% EndTime: 2019-12-29 18:25:49
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
	% StartTime: 2019-12-29 18:25:49
	% EndTime: 2019-12-29 18:25:49
	% DurationCPUTime: 0.13s
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:25:50
	% EndTime: 2019-12-29 18:25:50
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t134 = sin(qJ(2));
	t136 = cos(qJ(2));
	t146 = r_i_i_C(3) + qJ(3);
	t148 = pkin(2) + r_i_i_C(1);
	t149 = t148 * t134 - t146 * t136;
	t150 = t149 * qJD(2) - t134 * qJD(3);
	t147 = pkin(6) + r_i_i_C(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(2) * t137;
	t141 = -t146 * t134 - t148 * t136;
	t139 = -pkin(1) + t141;
	t1 = [t150 * t135 + (-t147 * t135 + t139 * t137) * qJD(1), (-t146 * t143 + t148 * t145) * t134 + (-t146 * t145 + (-t148 * qJD(2) + qJD(3)) * t137) * t136, -t134 * t145 + t136 * t143, 0, 0; -t150 * t137 + (t139 * t135 + t147 * t137) * qJD(1), -t149 * t144 + (t141 * qJD(2) + qJD(3) * t136) * t135, t135 * qJD(2) * t136 + t134 * t144, 0, 0; 0, -t150, qJD(2) * t134, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:25:49
	% EndTime: 2019-12-29 18:25:49
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (77->42), mult. (230->67), div. (0->0), fcn. (185->6), ass. (0->27)
	t36 = sin(qJ(2));
	t50 = t36 * qJD(3);
	t38 = cos(qJ(2));
	t55 = pkin(2) + pkin(3);
	t56 = -t38 * qJ(3) + t55 * t36;
	t57 = t56 * qJD(2) - t50;
	t37 = sin(qJ(1));
	t53 = qJD(1) * t37;
	t39 = cos(qJ(1));
	t52 = qJD(1) * t39;
	t51 = qJD(2) * t39;
	t49 = pkin(6) - r_i_i_C(3) - qJ(4);
	t34 = sin(pkin(8));
	t35 = cos(pkin(8));
	t48 = t34 * t38 - t35 * t36;
	t47 = t34 * t36 + t35 * t38;
	t46 = -qJ(3) * t36 - t55 * t38;
	t44 = t47 * t37;
	t43 = t47 * t39;
	t42 = qJD(1) * t48;
	t41 = t48 * qJD(2);
	t40 = -pkin(1) + t46;
	t33 = qJD(1) * t43 + t37 * t41;
	t32 = -qJD(2) * t44 + t39 * t42;
	t31 = -qJD(1) * t44 + t39 * t41;
	t30 = qJD(2) * t43 + t37 * t42;
	t1 = [-t33 * r_i_i_C(1) + t32 * r_i_i_C(2) - t39 * qJD(4) + t57 * t37 + (-t49 * t37 + t40 * t39) * qJD(1), -t30 * r_i_i_C(1) + t31 * r_i_i_C(2) + (-qJ(3) * t51 + t55 * t53) * t36 + (-qJ(3) * t53 + (-t55 * qJD(2) + qJD(3)) * t39) * t38, -t36 * t53 + t38 * t51, -t52, 0; t31 * r_i_i_C(1) + t30 * r_i_i_C(2) - t37 * qJD(4) - t57 * t39 + (t40 * t37 + t49 * t39) * qJD(1), t32 * r_i_i_C(1) + t33 * r_i_i_C(2) - t56 * t52 + (t46 * qJD(2) + qJD(3) * t38) * t37, t37 * qJD(2) * t38 + t36 * t52, -t53, 0; 0, t50 + ((t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + qJ(3)) * t38 + (-t35 * r_i_i_C(1) + t34 * r_i_i_C(2) - t55) * t36) * qJD(2), qJD(2) * t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:25:49
	% EndTime: 2019-12-29 18:25:49
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (234->45), mult. (424->68), div. (0->0), fcn. (364->8), ass. (0->35)
	t67 = pkin(8) + qJ(5);
	t65 = sin(t67);
	t66 = cos(t67);
	t70 = sin(qJ(2));
	t72 = cos(qJ(2));
	t81 = t65 * t72 - t66 * t70;
	t85 = pkin(4) * sin(pkin(8)) + qJ(3);
	t99 = pkin(2) + cos(pkin(8)) * pkin(4) + pkin(3);
	t100 = t99 * t70 - t85 * t72;
	t103 = t100 * qJD(2) - t70 * qJD(3);
	t102 = qJD(2) - qJD(5);
	t80 = t65 * t70 + t66 * t72;
	t59 = t102 * t80;
	t94 = qJD(2) * t70;
	t101 = t81 * qJD(5) + t66 * t94;
	t71 = sin(qJ(1));
	t96 = qJD(1) * t71;
	t73 = cos(qJ(1));
	t95 = qJD(1) * t73;
	t93 = qJD(2) * t72;
	t92 = qJD(2) * t73;
	t90 = pkin(6) - r_i_i_C(3) - pkin(7) - qJ(4);
	t86 = t72 * t92;
	t79 = qJD(1) * t81;
	t54 = t59 * t73 + t71 * t79;
	t78 = qJD(1) * t80;
	t55 = t101 * t73 - t65 * t86 + t71 * t78;
	t84 = t54 * r_i_i_C(1) + t55 * r_i_i_C(2);
	t56 = -t71 * t59 + t73 * t79;
	t57 = t73 * t78 + (t65 * t93 - t101) * t71;
	t83 = -t56 * r_i_i_C(1) - t57 * r_i_i_C(2);
	t82 = -t102 * t81 * r_i_i_C(1) - t59 * r_i_i_C(2);
	t77 = -t85 * t70 - t99 * t72;
	t75 = -pkin(1) + t77;
	t1 = [-t57 * r_i_i_C(1) + t56 * r_i_i_C(2) - t73 * qJD(4) + t103 * t71 + (-t90 * t71 + t75 * t73) * qJD(1), (-t85 * t92 + t99 * t96) * t70 + (-t85 * t96 + (-t99 * qJD(2) + qJD(3)) * t73) * t72 - t84, -t70 * t96 + t86, -t95, t84; -t55 * r_i_i_C(1) + t54 * r_i_i_C(2) - t71 * qJD(4) - t103 * t73 + (t75 * t71 + t90 * t73) * qJD(1), -t100 * t95 + (t77 * qJD(2) + qJD(3) * t72) * t71 - t83, t70 * t95 + t71 * t93, -t96, t83; 0, -t103 - t82, t94, 0, t82;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end