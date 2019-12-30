% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:24
	% EndTime: 2019-12-29 12:45:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:25
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:24
	% EndTime: 2019-12-29 12:45:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(6)) + r_i_i_C(2) * sin(pkin(6)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:25
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->12), mult. (48->17), div. (0->0), fcn. (34->4), ass. (0->9)
	t98 = r_i_i_C(2) + qJ(2);
	t92 = sin(qJ(1));
	t97 = qJD(1) * t92;
	t93 = cos(qJ(1));
	t96 = qJD(1) * t93;
	t90 = sin(pkin(6));
	t95 = qJD(3) * t90;
	t94 = -pkin(1) + (-pkin(2) - r_i_i_C(1)) * cos(pkin(6)) + (-r_i_i_C(3) - qJ(3)) * t90;
	t1 = [-t92 * t95 + t93 * qJD(2) + (-t98 * t92 + t94 * t93) * qJD(1), t96, -t90 * t97, 0; t93 * t95 + t92 * qJD(2) + (t94 * t92 + t98 * t93) * qJD(1), t97, t90 * t96, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:25
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (50->24), mult. (160->44), div. (0->0), fcn. (136->6), ass. (0->22)
	t36 = sin(qJ(1));
	t49 = qJD(1) * t36;
	t38 = cos(qJ(1));
	t48 = qJD(1) * t38;
	t33 = sin(pkin(6));
	t47 = qJD(3) * t33;
	t46 = pkin(5) + r_i_i_C(3) - qJ(2);
	t34 = cos(pkin(6));
	t35 = sin(qJ(4));
	t37 = cos(qJ(4));
	t45 = -t33 * t37 + t34 * t35;
	t44 = t33 * t35 + t34 * t37;
	t43 = t44 * t38;
	t42 = t44 * t36;
	t41 = qJD(1) * t45;
	t40 = t45 * qJD(4);
	t39 = -qJ(3) * t33 - pkin(1) + (-pkin(2) - pkin(3)) * t34;
	t32 = -qJD(1) * t43 + t36 * t40;
	t31 = qJD(4) * t42 + t38 * t41;
	t30 = qJD(1) * t42 + t38 * t40;
	t29 = -qJD(4) * t43 + t36 * t41;
	t1 = [-t36 * t47 + t32 * r_i_i_C(1) + t31 * r_i_i_C(2) + t38 * qJD(2) + (t46 * t36 + t39 * t38) * qJD(1), t48, -t33 * t49, r_i_i_C(1) * t29 + r_i_i_C(2) * t30; t38 * t47 - t30 * r_i_i_C(1) + t29 * r_i_i_C(2) + t36 * qJD(2) + (t39 * t36 - t46 * t38) * qJD(1), t49, t33 * t48, -r_i_i_C(1) * t31 + r_i_i_C(2) * t32; 0, 0, 0, (t45 * r_i_i_C(1) + t44 * r_i_i_C(2)) * qJD(4);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end