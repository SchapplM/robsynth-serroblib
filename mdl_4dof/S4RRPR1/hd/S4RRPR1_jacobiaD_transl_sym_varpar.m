% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR1
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:50:24
	% EndTime: 2019-10-09 20:50:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:50:24
	% EndTime: 2019-10-09 20:50:24
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
	% StartTime: 2019-10-09 20:50:24
	% EndTime: 2019-10-09 20:50:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:50:24
	% EndTime: 2019-10-09 20:50:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t49 = pkin(1) * qJD(1);
	t46 = qJ(1) + qJ(2);
	t42 = pkin(7) + t46;
	t40 = sin(t42);
	t41 = cos(t42);
	t45 = qJD(1) + qJD(2);
	t48 = (-pkin(2) * cos(t46) - r_i_i_C(1) * t41 + t40 * r_i_i_C(2)) * t45;
	t47 = (-pkin(2) * sin(t46) - r_i_i_C(1) * t40 - r_i_i_C(2) * t41) * t45;
	t1 = [-cos(qJ(1)) * t49 + t48, t48, 0, 0; -sin(qJ(1)) * t49 + t47, t47, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:50:24
	% EndTime: 2019-10-09 20:50:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (96->13), mult. (44->15), div. (0->0), fcn. (22->8), ass. (0->13)
	t60 = pkin(1) * qJD(1);
	t55 = qJ(1) + qJ(2);
	t54 = qJD(1) + qJD(2);
	t51 = pkin(7) + t55;
	t49 = qJ(4) + t51;
	t45 = sin(t49);
	t46 = cos(t49);
	t50 = qJD(4) + t54;
	t59 = (-r_i_i_C(1) * t46 + r_i_i_C(2) * t45) * t50;
	t58 = (-r_i_i_C(1) * t45 - r_i_i_C(2) * t46) * t50;
	t57 = (-pkin(2) * cos(t55) - pkin(3) * cos(t51)) * t54 + t59;
	t56 = (-pkin(2) * sin(t55) - pkin(3) * sin(t51)) * t54 + t58;
	t1 = [-cos(qJ(1)) * t60 + t57, t57, 0, t59; -sin(qJ(1)) * t60 + t56, t56, 0, t58; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end