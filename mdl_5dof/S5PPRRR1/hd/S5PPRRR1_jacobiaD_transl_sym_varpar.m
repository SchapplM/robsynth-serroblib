% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->5)
	t14 = pkin(9) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t17 = qJD(3) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, 0, cos(pkin(8)) * t17, 0, 0; 0, 0, sin(pkin(8)) * t17, 0, 0; 0, 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (48->10), mult. (40->14), div. (0->0), fcn. (25->6), ass. (0->15)
	t37 = pkin(3) * qJD(3);
	t29 = pkin(9) + qJ(3);
	t28 = qJ(4) + t29;
	t26 = cos(t28);
	t30 = qJD(3) + qJD(4);
	t36 = r_i_i_C(1) * t26 * t30;
	t25 = sin(t28);
	t35 = r_i_i_C(2) * t25 * t30;
	t34 = (-r_i_i_C(1) * t25 - r_i_i_C(2) * t26) * t30;
	t33 = -cos(t29) * t37 - t36;
	t32 = cos(pkin(8));
	t31 = sin(pkin(8));
	t24 = t32 * t35;
	t23 = t31 * t35;
	t1 = [0, 0, t33 * t32 + t24, -t32 * t36 + t24, 0; 0, 0, t33 * t31 + t23, -t31 * t36 + t23, 0; 0, 0, -sin(t29) * t37 + t34, t34, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:12
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (199->30), mult. (192->54), div. (0->0), fcn. (139->8), ass. (0->30)
	t191 = pkin(9) + qJ(3);
	t190 = qJ(4) + t191;
	t187 = sin(t190);
	t188 = cos(t190);
	t196 = cos(qJ(5));
	t209 = qJD(5) * t196;
	t192 = qJD(3) + qJD(4);
	t195 = sin(qJ(5));
	t216 = t192 * t195;
	t222 = t187 * t209 + t188 * t216;
	t221 = pkin(7) + r_i_i_C(3);
	t210 = qJD(5) * t195;
	t205 = t187 * t210;
	t220 = r_i_i_C(1) * t205 + t222 * r_i_i_C(2);
	t219 = pkin(3) * qJD(3);
	t218 = t187 * t192;
	t215 = t192 * t196;
	t193 = sin(pkin(8));
	t214 = t193 * t195;
	t213 = t193 * t196;
	t194 = cos(pkin(8));
	t212 = t194 * t195;
	t211 = t194 * t196;
	t207 = t220 * t193;
	t206 = t220 * t194;
	t200 = (r_i_i_C(1) * t195 + r_i_i_C(2) * t196) * t218;
	t199 = ((-r_i_i_C(1) * t196 - pkin(4)) * t188 - t221 * t187) * t192;
	t198 = -cos(t191) * t219 + t199;
	t197 = -pkin(4) * t218 + (-t187 * t215 - t188 * t210) * r_i_i_C(1) + t221 * t188 * t192 + (t187 * t216 - t188 * t209) * r_i_i_C(2);
	t1 = [0, 0, t198 * t194 + t206, t194 * t199 + t206, t194 * t200 + ((-t188 * t211 - t214) * r_i_i_C(1) + (t188 * t212 - t213) * r_i_i_C(2)) * qJD(5); 0, 0, t198 * t193 + t207, t193 * t199 + t207, t193 * t200 + ((-t188 * t213 + t212) * r_i_i_C(1) + (t188 * t214 + t211) * r_i_i_C(2)) * qJD(5); 0, 0, -sin(t191) * t219 + t197, t197, (-t188 * t215 + t205) * r_i_i_C(2) - t222 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end