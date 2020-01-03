% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
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
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (36->15), mult. (88->23), div. (0->0), fcn. (74->4), ass. (0->14)
	t87 = -pkin(1) - pkin(2);
	t78 = sin(qJ(1));
	t86 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t85 = qJD(1) * t80;
	t84 = qJD(3) * t78;
	t83 = qJD(3) * t80;
	t77 = sin(qJ(3));
	t79 = cos(qJ(3));
	t71 = -t77 * t84 - t79 * t83 + (t77 * t78 + t79 * t80) * qJD(1);
	t72 = (t84 - t86) * t79 + (-t83 + t85) * t77;
	t82 = -t72 * r_i_i_C(1) - t71 * r_i_i_C(2);
	t81 = t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
	t1 = [t80 * qJD(2) + (-qJ(2) * t78 + t87 * t80) * qJD(1) - t81, t85, t81, 0, 0; t78 * qJD(2) + (qJ(2) * t80 + t87 * t78) * qJD(1) - t82, t86, t82, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (84->25), mult. (132->35), div. (0->0), fcn. (104->6), ass. (0->22)
	t105 = (qJD(1) - qJD(3)) * pkin(3);
	t91 = cos(qJ(3));
	t103 = -t91 * pkin(3) - pkin(1) - pkin(2);
	t102 = pkin(3) * qJD(3);
	t90 = sin(qJ(1));
	t101 = qJD(1) * t90;
	t92 = cos(qJ(1));
	t100 = qJD(1) * t92;
	t99 = qJD(3) * t90;
	t98 = qJD(3) * t92;
	t89 = sin(qJ(3));
	t97 = pkin(3) * t89 + qJ(2);
	t88 = qJ(3) + pkin(8);
	t86 = sin(t88);
	t87 = cos(t88);
	t79 = -t86 * t99 - t87 * t98 + (t86 * t90 + t87 * t92) * qJD(1);
	t80 = (-t101 + t99) * t87 + (t100 - t98) * t86;
	t96 = -t80 * r_i_i_C(1) - t79 * r_i_i_C(2);
	t95 = t79 * r_i_i_C(1) - t80 * r_i_i_C(2);
	t94 = t89 * t92 - t90 * t91;
	t93 = t89 * t90 + t91 * t92;
	t1 = [t92 * qJD(2) + t93 * t102 + (t103 * t92 - t97 * t90) * qJD(1) - t95, t100, t93 * t105 + t95, 0, 0; t90 * qJD(2) - t94 * t102 + (t103 * t90 + t97 * t92) * qJD(1) - t96, t101, -t94 * t105 + t96, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (259->37), mult. (352->59), div. (0->0), fcn. (306->8), ass. (0->29)
	t97 = (qJD(1) - qJD(3)) * pkin(3);
	t67 = sin(qJ(1));
	t70 = cos(qJ(1));
	t85 = qJ(3) + pkin(8);
	t82 = sin(t85);
	t83 = cos(t85);
	t71 = t67 * t82 + t70 * t83;
	t74 = qJD(3) * t82;
	t75 = qJD(3) * t83;
	t53 = t71 * qJD(1) - t67 * t74 - t70 * t75;
	t56 = -t67 * t83 + t70 * t82;
	t54 = t56 * qJD(1) + t67 * t75 - t70 * t74;
	t65 = sin(qJ(5));
	t68 = cos(qJ(5));
	t76 = r_i_i_C(1) * t68 - r_i_i_C(2) * t65 + pkin(4);
	t90 = pkin(7) + r_i_i_C(3);
	t94 = (r_i_i_C(1) * t65 + r_i_i_C(2) * t68) * qJD(5);
	t96 = -t76 * t53 - t90 * t54 - t56 * t94;
	t95 = -t90 * t53 + t76 * t54 - t71 * t94;
	t69 = cos(qJ(3));
	t89 = -t69 * pkin(3) - pkin(1) - pkin(2);
	t88 = pkin(3) * qJD(3);
	t87 = qJD(5) * t65;
	t86 = qJD(5) * t68;
	t66 = sin(qJ(3));
	t84 = pkin(3) * t66 + qJ(2);
	t80 = t66 * t70 - t67 * t69;
	t79 = t66 * t67 + t69 * t70;
	t1 = [t70 * qJD(2) + t79 * t88 + (-t84 * t67 + t89 * t70) * qJD(1) + t96, qJD(1) * t70, t79 * t97 - t96, 0, (-t54 * t68 + t71 * t87) * r_i_i_C(2) + (-t54 * t65 - t71 * t86) * r_i_i_C(1); t67 * qJD(2) - t80 * t88 + (t89 * t67 + t84 * t70) * qJD(1) + t95, qJD(1) * t67, -t80 * t97 - t95, 0, (-t53 * t68 - t56 * t87) * r_i_i_C(2) + (-t53 * t65 + t56 * t86) * r_i_i_C(1); 0, 0, 0, 0, t94;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end