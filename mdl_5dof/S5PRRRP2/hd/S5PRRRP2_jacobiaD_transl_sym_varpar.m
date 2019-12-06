% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:42:26
	% EndTime: 2019-12-05 16:42:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:42:26
	% EndTime: 2019-12-05 16:42:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:42:26
	% EndTime: 2019-12-05 16:42:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:42:26
	% EndTime: 2019-12-05 16:42:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(8) + qJ(2);
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
	% StartTime: 2019-12-05 16:42:26
	% EndTime: 2019-12-05 16:42:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->20), mult. (104->33), div. (0->0), fcn. (64->6), ass. (0->20)
	t63 = r_i_i_C(3) + pkin(7);
	t43 = pkin(8) + qJ(2);
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
	% StartTime: 2019-12-05 16:42:27
	% EndTime: 2019-12-05 16:42:27
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (248->25), mult. (202->34), div. (0->0), fcn. (136->6), ass. (0->27)
	t174 = pkin(8) + qJ(2);
	t173 = qJ(3) + t174;
	t172 = cos(t173);
	t175 = qJD(2) + qJD(3);
	t192 = t172 * t175;
	t196 = pkin(4) + r_i_i_C(1);
	t197 = pkin(7) + r_i_i_C(2);
	t195 = r_i_i_C(3) + qJ(5);
	t177 = cos(qJ(4));
	t184 = t195 * t177;
	t176 = sin(qJ(4));
	t187 = t196 * t176;
	t182 = -t187 + t184;
	t194 = pkin(2) * qJD(2);
	t171 = sin(t173);
	t193 = t171 * t175;
	t191 = t175 * t176;
	t190 = qJD(4) * t176;
	t189 = qJD(4) * t177;
	t188 = t176 * qJD(5);
	t185 = t172 * t189;
	t183 = -t195 * t176 - t196 * t177;
	t181 = -pkin(3) + t183;
	t180 = t183 * qJD(4) + qJD(5) * t177;
	t179 = t181 * t193 + t197 * t192 + t195 * t185 + (-qJD(4) * t187 + t188) * t172;
	t178 = t181 * t192 + (-qJD(4) * t184 - t197 * t175 + t196 * t190 - t188) * t171;
	t1 = [0, -cos(t174) * t194 + t178, t178, t180 * t172 - t182 * t193, -t171 * t191 + t185; 0, -sin(t174) * t194 + t179, t179, t180 * t171 + t182 * t192, t171 * t189 + t172 * t191; 0, 0, 0, t182 * qJD(4) + t188, t190;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end