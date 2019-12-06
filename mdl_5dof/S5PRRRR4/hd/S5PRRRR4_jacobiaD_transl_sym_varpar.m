% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(9) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(9) + qJ(2);
	t40 = qJ(3) + t41;
	t38 = sin(t40);
	t39 = cos(t40);
	t42 = qJD(2) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t42;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t42;
	t1 = [0, -cos(t41) * t45 + t44, t44, 0, 0; 0, -sin(t41) * t45 + t43, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->20), mult. (104->33), div. (0->0), fcn. (64->6), ass. (0->20)
	t63 = r_i_i_C(3) + pkin(7);
	t43 = pkin(9) + qJ(2);
	t42 = qJ(3) + t43;
	t40 = sin(t42);
	t41 = cos(t42);
	t46 = cos(qJ(4));
	t55 = qJD(4) * t46;
	t44 = qJD(2) + qJD(3);
	t45 = sin(qJ(4));
	t59 = t44 * t45;
	t62 = t40 * t59 - t41 * t55;
	t61 = t40 * t55 + t41 * t59;
	t58 = t44 * t46;
	t57 = pkin(2) * qJD(2);
	t56 = qJD(4) * t45;
	t52 = t40 * t56;
	t49 = t40 * t58 + t41 * t56;
	t48 = r_i_i_C(1) * t52 + ((-r_i_i_C(1) * t46 - pkin(3)) * t41 - t63 * t40) * t44 + t61 * r_i_i_C(2);
	t47 = -r_i_i_C(1) * t49 + t62 * r_i_i_C(2) + (-pkin(3) * t40 + t41 * t63) * t44;
	t1 = [0, -cos(t43) * t57 + t48, t48, t62 * r_i_i_C(1) + t49 * r_i_i_C(2), 0; 0, -sin(t43) * t57 + t47, t47, (-t41 * t58 + t52) * r_i_i_C(2) - t61 * r_i_i_C(1), 0; 0, 0, 0, (-r_i_i_C(1) * t45 - r_i_i_C(2) * t46) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:08:31
	% EndTime: 2019-12-05 17:08:31
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (249->33), mult. (162->45), div. (0->0), fcn. (103->8), ass. (0->36)
	t61 = pkin(9) + qJ(2);
	t58 = qJ(3) + t61;
	t56 = cos(t58);
	t63 = qJD(2) + qJD(3);
	t88 = t56 * t63;
	t55 = sin(t58);
	t64 = qJ(4) + qJ(5);
	t60 = cos(t64);
	t62 = qJD(4) + qJD(5);
	t85 = t60 * t62;
	t59 = sin(t64);
	t86 = t59 * t63;
	t92 = t55 * t85 + t56 * t86;
	t65 = sin(qJ(4));
	t91 = pkin(4) * t65;
	t90 = r_i_i_C(2) * t60;
	t89 = t55 * t63;
	t87 = t59 * t62;
	t84 = t63 * (-pkin(8) - pkin(7));
	t83 = pkin(2) * qJD(2);
	t66 = cos(qJ(4));
	t82 = qJD(4) * t66;
	t81 = r_i_i_C(1) * t85;
	t80 = t63 * t90;
	t79 = t55 * t87;
	t78 = t55 * t86;
	t75 = qJD(4) * t91;
	t74 = -pkin(4) * t66 - r_i_i_C(1) * t60 - pkin(3);
	t73 = -r_i_i_C(1) * t59 - t90;
	t72 = t73 * t62;
	t71 = r_i_i_C(1) * t78 + t55 * t80 + (r_i_i_C(2) * t87 - t81) * t56;
	t70 = t72 - t75;
	t69 = t74 * t88 + r_i_i_C(1) * t79 + (-r_i_i_C(3) * t63 + t75 + t84) * t55 + t92 * r_i_i_C(2);
	t68 = r_i_i_C(2) * t78 + r_i_i_C(3) * t88 + t74 * t89 + (t70 - t84) * t56;
	t44 = r_i_i_C(2) * t79;
	t1 = [0, -cos(t61) * t83 + t69, t69, (-t56 * t82 + t65 * t89) * pkin(4) + t71, t71; 0, -sin(t61) * t83 + t68, t68, t44 + (-pkin(4) * t82 - t81) * t55 + (t73 - t91) * t88, -t92 * r_i_i_C(1) - t56 * t80 + t44; 0, 0, 0, t70, t72;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end