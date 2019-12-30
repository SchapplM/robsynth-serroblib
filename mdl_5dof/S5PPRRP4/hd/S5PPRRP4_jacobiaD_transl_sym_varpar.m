% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:23
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:15
	% EndTime: 2019-12-29 15:23:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:09
	% EndTime: 2019-12-29 15:23:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:21
	% EndTime: 2019-12-29 15:23:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:21
	% EndTime: 2019-12-29 15:23:21
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
	% StartTime: 2019-12-29 15:23:16
	% EndTime: 2019-12-29 15:23:16
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (41->16), mult. (122->31), div. (0->0), fcn. (106->6), ass. (0->16)
	t45 = cos(pkin(7));
	t47 = sin(qJ(3));
	t55 = sin(pkin(7));
	t57 = cos(qJ(3));
	t41 = -t45 * t47 + t55 * t57;
	t46 = sin(qJ(4));
	t48 = cos(qJ(4));
	t49 = (r_i_i_C(1) * t46 + r_i_i_C(2) * t48) * qJD(4);
	t58 = -pkin(6) - r_i_i_C(3);
	t54 = qJD(4) * t46;
	t53 = qJD(4) * t48;
	t50 = r_i_i_C(1) * t48 - r_i_i_C(2) * t46 + pkin(3);
	t40 = -t45 * t57 - t55 * t47;
	t39 = t40 * qJD(3);
	t38 = t41 * qJD(3);
	t1 = [0, 0, -t58 * t38 + t50 * t39 - t41 * t49, (-t38 * t48 - t40 * t54) * r_i_i_C(2) + (-t38 * t46 + t40 * t53) * r_i_i_C(1), 0; 0, 0, -t50 * t38 - t58 * t39 - t40 * t49, (-t39 * t48 + t41 * t54) * r_i_i_C(2) + (-t39 * t46 - t41 * t53) * r_i_i_C(1), 0; 0, 0, 0, t49, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:15
	% EndTime: 2019-12-29 15:23:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (64->17), mult. (172->24), div. (0->0), fcn. (155->6), ass. (0->18)
	t53 = sin(qJ(4));
	t55 = cos(qJ(4));
	t67 = pkin(4) + r_i_i_C(1);
	t68 = -r_i_i_C(2) * t53 + t67 * t55;
	t51 = cos(pkin(7));
	t54 = sin(qJ(3));
	t61 = sin(pkin(7));
	t63 = cos(qJ(3));
	t45 = t51 * t54 - t61 * t63;
	t58 = r_i_i_C(2) * t55 + t67 * t53;
	t56 = t58 * qJD(4);
	t64 = -r_i_i_C(3) - qJ(5) - pkin(6);
	t59 = pkin(3) + t68;
	t57 = qJD(4) * t68;
	t44 = -t51 * t63 - t61 * t54;
	t43 = t44 * qJD(3);
	t42 = t45 * qJD(3);
	t1 = [0, 0, -t44 * qJD(5) + t64 * t42 + t59 * t43 + t45 * t56, t58 * t42 + t44 * t57, -t43; 0, 0, -t45 * qJD(5) + t59 * t42 - t64 * t43 - t44 * t56, -t58 * t43 + t45 * t57, -t42; 0, 0, 0, t56, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end