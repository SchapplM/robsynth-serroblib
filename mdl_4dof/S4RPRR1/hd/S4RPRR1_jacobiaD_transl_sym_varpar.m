% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:46:09
	% EndTime: 2019-10-09 20:46:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:46:09
	% EndTime: 2019-10-09 20:46:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:46:09
	% EndTime: 2019-10-09 20:46:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(7);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:46:09
	% EndTime: 2019-10-09 20:46:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t42 = qJ(1) + pkin(7);
	t40 = qJ(3) + t42;
	t38 = sin(t40);
	t39 = cos(t40);
	t41 = qJD(1) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t41;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t41;
	t1 = [(-cos(t42) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t44, 0, t44, 0; t43 + (-sin(t42) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t43, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:46:09
	% EndTime: 2019-10-09 20:46:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (88->13), mult. (40->15), div. (0->0), fcn. (20->8), ass. (0->13)
	t50 = qJD(1) + qJD(3);
	t56 = pkin(3) * t50;
	t51 = qJ(1) + pkin(7);
	t49 = qJ(3) + t51;
	t47 = qJ(4) + t49;
	t43 = sin(t47);
	t44 = cos(t47);
	t48 = qJD(4) + t50;
	t55 = (-r_i_i_C(1) * t44 + r_i_i_C(2) * t43) * t48;
	t54 = (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * t48;
	t53 = -cos(t49) * t56 + t55;
	t52 = -sin(t49) * t56 + t54;
	t1 = [(-cos(t51) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t53, 0, t53, t55; (-sin(t51) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t52, 0, t52, t54; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end