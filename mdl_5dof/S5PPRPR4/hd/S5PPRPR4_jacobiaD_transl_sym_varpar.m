% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (6->4), mult. (20->10), div. (0->0), fcn. (16->4), ass. (0->7)
	t49 = cos(qJ(3));
	t48 = sin(qJ(3));
	t47 = cos(pkin(7));
	t46 = sin(pkin(7));
	t45 = (t46 * t48 + t47 * t49) * qJD(3);
	t44 = (-t46 * t49 + t47 * t48) * qJD(3);
	t1 = [0, 0, -t45 * r_i_i_C(1) + t44 * r_i_i_C(2), 0, 0; 0, 0, t44 * r_i_i_C(1) + t45 * r_i_i_C(2), 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:54
	% EndTime: 2019-12-29 15:17:54
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (24->10), mult. (68->14), div. (0->0), fcn. (60->6), ass. (0->11)
	t26 = sin(pkin(7));
	t28 = cos(pkin(7));
	t29 = sin(qJ(3));
	t30 = cos(qJ(3));
	t36 = -t26 * t30 + t28 * t29;
	t33 = -r_i_i_C(3) - qJ(4);
	t32 = t26 * t29 + t28 * t30;
	t31 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(3);
	t23 = t32 * qJD(3);
	t22 = t36 * qJD(3);
	t1 = [0, 0, t32 * qJD(4) + t33 * t22 - t31 * t23, t23, 0; 0, 0, -t36 * qJD(4) + t31 * t22 + t33 * t23, -t22, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:54
	% EndTime: 2019-12-29 15:17:54
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (69->22), mult. (138->34), div. (0->0), fcn. (124->7), ass. (0->17)
	t54 = cos(pkin(7));
	t56 = sin(qJ(3));
	t63 = sin(pkin(7));
	t65 = cos(qJ(3));
	t45 = t54 * t56 - t63 * t65;
	t53 = pkin(8) + qJ(5);
	t51 = sin(t53);
	t52 = cos(t53);
	t57 = (r_i_i_C(1) * t51 + r_i_i_C(2) * t52) * qJD(5);
	t66 = -r_i_i_C(3) - pkin(6) - qJ(4);
	t62 = qJD(5) * t51;
	t61 = qJD(5) * t52;
	t58 = r_i_i_C(1) * t52 - r_i_i_C(2) * t51 + cos(pkin(8)) * pkin(4) + pkin(3);
	t44 = -t54 * t65 - t56 * t63;
	t43 = t44 * qJD(3);
	t42 = t45 * qJD(3);
	t1 = [0, 0, -t44 * qJD(4) + t42 * t66 + t58 * t43 + t45 * t57, -t43, (t42 * t52 - t44 * t62) * r_i_i_C(2) + (t42 * t51 + t44 * t61) * r_i_i_C(1); 0, 0, -t45 * qJD(4) + t58 * t42 - t43 * t66 - t44 * t57, -t42, (-t43 * t52 - t45 * t62) * r_i_i_C(2) + (-t43 * t51 + t45 * t61) * r_i_i_C(1); 0, 0, 0, 0, t57;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end