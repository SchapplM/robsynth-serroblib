% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:33:23
	% EndTime: 2019-10-09 20:33:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:33:23
	% EndTime: 2019-10-09 20:33:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:33:23
	% EndTime: 2019-10-09 20:33:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(6) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:33:23
	% EndTime: 2019-10-09 20:33:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(6) + qJ(2);
	t40 = qJ(3) + t41;
	t38 = sin(t40);
	t39 = cos(t40);
	t42 = qJD(2) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t42;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t42;
	t1 = [0, -cos(t41) * t45 + t44, t44, 0; 0, -sin(t41) * t45 + t43, t43, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:33:23
	% EndTime: 2019-10-09 20:33:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->11), mult. (42->12), div. (0->0), fcn. (24->4), ass. (0->13)
	t33 = qJ(4) + r_i_i_C(3);
	t32 = -pkin(3) - r_i_i_C(1);
	t25 = pkin(6) + qJ(2);
	t24 = qJ(3) + t25;
	t22 = sin(t24);
	t26 = qJD(2) + qJD(3);
	t31 = t26 * t22;
	t23 = cos(t24);
	t30 = t26 * t23;
	t29 = pkin(2) * qJD(2);
	t28 = t22 * qJD(4) + t33 * t30 + t32 * t31;
	t27 = t23 * qJD(4) + (-t33 * t22 + t32 * t23) * t26;
	t1 = [0, -cos(t25) * t29 + t27, t27, t30; 0, -sin(t25) * t29 + t28, t28, t31; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end