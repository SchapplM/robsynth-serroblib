% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(2));
	t26 = sin(qJ(2));
	t1 = [0, (-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(2) * qJD(2);
	t40 = qJ(2) + qJ(3);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(2) + qJD(3);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [0, -cos(qJ(2)) * t43 + t42, t42, 0, 0; 0, -sin(qJ(2)) * t43 + t41, t41, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(2) + qJD(3);
	t55 = pkin(3) * t48;
	t54 = pkin(2) * qJD(2);
	t49 = qJ(2) + qJ(3);
	t47 = qJ(4) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t44 = qJD(4) + t48;
	t53 = (-r_i_i_C(1) * t43 + r_i_i_C(2) * t42) * t44;
	t52 = (-r_i_i_C(1) * t42 - r_i_i_C(2) * t43) * t44;
	t51 = -cos(t49) * t55 + t53;
	t50 = -sin(t49) * t55 + t52;
	t1 = [0, -cos(qJ(2)) * t54 + t51, t51, t53, 0; 0, -sin(qJ(2)) * t54 + t50, t50, t52, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:56:14
	% EndTime: 2019-10-09 20:56:14
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (195->21), mult. (136->35), div. (0->0), fcn. (84->8), ass. (0->24)
	t71 = pkin(6) + r_i_i_C(3);
	t49 = qJ(2) + qJ(3);
	t47 = qJ(4) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t51 = cos(qJ(5));
	t62 = qJD(5) * t51;
	t48 = qJD(2) + qJD(3);
	t44 = qJD(4) + t48;
	t50 = sin(qJ(5));
	t66 = t44 * t50;
	t70 = t42 * t66 - t43 * t62;
	t69 = t42 * t62 + t43 * t66;
	t68 = pkin(3) * t48;
	t65 = t44 * t51;
	t64 = pkin(2) * qJD(2);
	t63 = qJD(5) * t50;
	t59 = t42 * t63;
	t56 = t42 * t65 + t43 * t63;
	t55 = r_i_i_C(1) * t59 + (-r_i_i_C(1) * t43 * t51 - t71 * t42) * t44 + t69 * r_i_i_C(2);
	t54 = -cos(t49) * t68 + t55;
	t53 = t71 * t43 * t44 - t56 * r_i_i_C(1) + t70 * r_i_i_C(2);
	t52 = -sin(t49) * t68 + t53;
	t1 = [0, -cos(qJ(2)) * t64 + t54, t54, t55, t70 * r_i_i_C(1) + t56 * r_i_i_C(2); 0, -sin(qJ(2)) * t64 + t52, t52, t53, (-t43 * t65 + t59) * r_i_i_C(2) - t69 * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t50 - r_i_i_C(2) * t51) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end