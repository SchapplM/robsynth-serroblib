% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:28:35
	% EndTime: 2019-12-31 18:28:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:28:35
	% EndTime: 2019-12-31 18:28:35
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
	% StartTime: 2019-12-31 18:28:35
	% EndTime: 2019-12-31 18:28:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:28:35
	% EndTime: 2019-12-31 18:28:35
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(6) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(8) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(8)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t28, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:28:36
	% EndTime: 2019-12-31 18:28:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (96->25), mult. (140->39), div. (0->0), fcn. (98->5), ass. (0->16)
	t142 = pkin(8) + qJ(3);
	t140 = sin(t142);
	t141 = cos(t142);
	t154 = r_i_i_C(3) + qJ(4);
	t156 = pkin(3) + r_i_i_C(1);
	t158 = (t156 * t140 - t154 * t141) * qJD(3) - t140 * qJD(4);
	t155 = r_i_i_C(2) + pkin(6) + qJ(2);
	t144 = sin(qJ(1));
	t153 = qJD(1) * t144;
	t145 = cos(qJ(1));
	t152 = qJD(1) * t145;
	t151 = qJD(3) * t141;
	t149 = qJD(3) * t154;
	t148 = -t156 * qJD(3) + qJD(4);
	t147 = -t154 * t140 - t156 * t141 - cos(pkin(8)) * pkin(2) - pkin(1);
	t1 = [t145 * qJD(2) + t158 * t144 + (-t155 * t144 + t147 * t145) * qJD(1), t152, (-t145 * t149 + t156 * t153) * t140 + (t148 * t145 - t154 * t153) * t141, -t140 * t153 + t145 * t151, 0; t144 * qJD(2) - t158 * t145 + (t147 * t144 + t155 * t145) * qJD(1), t153, (-t144 * t149 - t156 * t152) * t140 + (t148 * t144 + t154 * t152) * t141, t140 * t152 + t144 * t151, 0; 0, 0, -t158, qJD(3) * t140, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:28:35
	% EndTime: 2019-12-31 18:28:36
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (253->43), mult. (382->70), div. (0->0), fcn. (331->7), ass. (0->34)
	t67 = pkin(8) + qJ(3);
	t65 = sin(t67);
	t66 = cos(t67);
	t69 = sin(qJ(5));
	t71 = cos(qJ(5));
	t103 = -t65 * t71 + t66 * t69;
	t101 = qJD(3) - qJD(5);
	t78 = t65 * t69 + t66 * t71;
	t59 = t101 * t78;
	t98 = pkin(3) + pkin(4);
	t102 = (-qJ(4) * t66 + t98 * t65) * qJD(3) - t65 * qJD(4);
	t92 = qJD(3) * t65;
	t100 = t103 * qJD(5) + t71 * t92;
	t70 = sin(qJ(1));
	t94 = qJD(1) * t70;
	t72 = cos(qJ(1));
	t93 = qJD(1) * t72;
	t91 = qJD(3) * t66;
	t89 = qJ(4) * qJD(1);
	t88 = qJ(4) * qJD(3);
	t87 = pkin(7) + r_i_i_C(3) - pkin(6) - qJ(2);
	t83 = t72 * t91;
	t76 = qJD(1) * t103;
	t54 = t59 * t72 + t70 * t76;
	t75 = qJD(1) * t78;
	t55 = t100 * t72 - t69 * t83 + t70 * t75;
	t82 = t54 * r_i_i_C(1) + t55 * r_i_i_C(2);
	t56 = -t70 * t59 + t72 * t76;
	t57 = t72 * t75 + (t69 * t91 - t100) * t70;
	t81 = -t56 * r_i_i_C(1) - t57 * r_i_i_C(2);
	t80 = -r_i_i_C(1) * t101 * t103 - t59 * r_i_i_C(2);
	t77 = -t98 * qJD(3) + qJD(4);
	t74 = -qJ(4) * t65 - t98 * t66 - cos(pkin(8)) * pkin(2) - pkin(1);
	t1 = [-t57 * r_i_i_C(1) + t56 * r_i_i_C(2) + t72 * qJD(2) + t102 * t70 + (t87 * t70 + t74 * t72) * qJD(1), t93, (-t72 * t88 + t98 * t94) * t65 + (-t70 * t89 + t77 * t72) * t66 - t82, -t65 * t94 + t83, t82; -t55 * r_i_i_C(1) + t54 * r_i_i_C(2) + t70 * qJD(2) - t102 * t72 + (t74 * t70 - t87 * t72) * qJD(1), t94, (-t70 * t88 - t98 * t93) * t65 + (t77 * t70 + t72 * t89) * t66 - t81, t65 * t93 + t70 * t91, t81; 0, 0, -t102 - t80, t92, t80;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end