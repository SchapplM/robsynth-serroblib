% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:44
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
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:44
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:44
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(8));
	t54 = sin(pkin(8));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t57 * t58) * qJD(1), qJD(1) * t57, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t56 * t58) * qJD(1), qJD(1) * t56, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (72->18), mult. (96->25), div. (0->0), fcn. (80->6), ass. (0->16)
	t98 = -pkin(1) - cos(pkin(8)) * pkin(3) - pkin(2);
	t89 = sin(qJ(1));
	t97 = qJD(1) * t89;
	t90 = cos(qJ(1));
	t96 = qJD(1) * t90;
	t95 = qJD(4) * t89;
	t94 = qJD(4) * t90;
	t93 = pkin(3) * sin(pkin(8)) + qJ(2);
	t87 = pkin(8) + qJ(4);
	t85 = sin(t87);
	t86 = cos(t87);
	t78 = -t85 * t95 - t86 * t94 + (t85 * t89 + t86 * t90) * qJD(1);
	t79 = (t95 - t97) * t86 + (-t94 + t96) * t85;
	t92 = -t79 * r_i_i_C(1) - t78 * r_i_i_C(2);
	t91 = t78 * r_i_i_C(1) - t79 * r_i_i_C(2);
	t1 = [t90 * qJD(2) + (-t93 * t89 + t98 * t90) * qJD(1) - t91, t96, 0, t91, 0; t89 * qJD(2) + (t98 * t89 + t93 * t90) * qJD(1) - t92, t97, 0, t92, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:18:44
	% EndTime: 2019-12-29 16:18:45
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (247->26), mult. (316->40), div. (0->0), fcn. (282->8), ass. (0->24)
	t94 = qJD(1) - qJD(4);
	t78 = pkin(8) + qJ(4);
	t76 = sin(t78);
	t77 = cos(t78);
	t81 = sin(qJ(1));
	t82 = cos(qJ(1));
	t55 = -t81 * t76 - t82 * t77;
	t53 = t94 * t55;
	t56 = t82 * t76 - t81 * t77;
	t54 = t94 * t56;
	t66 = sin(qJ(5));
	t67 = cos(qJ(5));
	t74 = r_i_i_C(1) * t67 - r_i_i_C(2) * t66 + pkin(4);
	t84 = pkin(7) + r_i_i_C(3);
	t88 = (r_i_i_C(1) * t66 + r_i_i_C(2) * t67) * qJD(5);
	t93 = -t84 * t53 - t74 * t54 - t55 * t88;
	t92 = -t74 * t53 + t84 * t54 + t56 * t88;
	t91 = qJ(2) + pkin(3) * sin(pkin(8));
	t90 = qJD(1) * t81;
	t89 = qJD(1) * t82;
	t85 = qJD(2) + (-pkin(1) - cos(pkin(8)) * pkin(3) - pkin(2)) * qJD(1);
	t80 = qJD(5) * t66;
	t79 = qJD(5) * t67;
	t1 = [t85 * t82 - t91 * t90 - t92, t89, 0, t92, (-t54 * t67 - t55 * t80) * r_i_i_C(2) + (-t54 * t66 + t55 * t79) * r_i_i_C(1); t85 * t81 + t91 * t89 - t93, t90, 0, t93, (t53 * t67 - t56 * t80) * r_i_i_C(2) + (t53 * t66 + t56 * t79) * r_i_i_C(1); 0, 0, 0, 0, t88;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end