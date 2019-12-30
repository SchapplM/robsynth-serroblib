% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(8)) * t14, 0, 0, 0; 0, sin(pkin(8)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (33->9), mult. (40->14), div. (0->0), fcn. (25->6), ass. (0->14)
	t35 = pkin(2) * qJD(2);
	t27 = qJ(2) + qJ(3);
	t25 = cos(t27);
	t26 = qJD(2) + qJD(3);
	t34 = r_i_i_C(1) * t25 * t26;
	t24 = sin(t27);
	t33 = r_i_i_C(2) * t24 * t26;
	t32 = (-r_i_i_C(1) * t24 - r_i_i_C(2) * t25) * t26;
	t31 = -cos(qJ(2)) * t35 - t34;
	t29 = cos(pkin(8));
	t28 = sin(pkin(8));
	t23 = t29 * t33;
	t22 = t28 * t33;
	t1 = [0, t31 * t29 + t23, -t29 * t34 + t23, 0, 0; 0, t31 * t28 + t22, -t28 * t34 + t22, 0, 0; 0, -sin(qJ(2)) * t35 + t32, t32, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:22
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (63->12), mult. (54->17), div. (0->0), fcn. (33->8), ass. (0->15)
	t36 = qJ(2) + qJ(3);
	t32 = pkin(9) + t36;
	t31 = cos(t32);
	t35 = qJD(2) + qJD(3);
	t40 = t35 * (-r_i_i_C(1) * t31 - pkin(3) * cos(t36));
	t43 = pkin(2) * qJD(2);
	t30 = sin(t32);
	t42 = r_i_i_C(2) * t30 * t35;
	t41 = -cos(qJ(2)) * t43 + t40;
	t39 = (-pkin(3) * sin(t36) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * t35;
	t38 = cos(pkin(8));
	t37 = sin(pkin(8));
	t29 = t38 * t42;
	t28 = t37 * t42;
	t1 = [0, t41 * t38 + t29, t38 * t40 + t29, 0, 0; 0, t41 * t37 + t28, t37 * t40 + t28, 0, 0; 0, -sin(qJ(2)) * t43 + t39, t39, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:23
	% EndTime: 2019-12-29 15:40:23
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (214->30), mult. (206->52), div. (0->0), fcn. (147->10), ass. (0->30)
	t198 = qJ(2) + qJ(3);
	t194 = pkin(9) + t198;
	t193 = cos(t194);
	t202 = cos(qJ(5));
	t192 = sin(t194);
	t218 = qJD(5) * t192;
	t197 = qJD(2) + qJD(3);
	t201 = sin(qJ(5));
	t223 = t197 * t201;
	t229 = t193 * t223 + t202 * t218;
	t228 = pkin(7) + r_i_i_C(3);
	t213 = t201 * t218;
	t227 = r_i_i_C(1) * t213 + t229 * r_i_i_C(2);
	t214 = -r_i_i_C(1) * t202 - pkin(4);
	t204 = t197 * (-t228 * t192 + t214 * t193 - pkin(3) * cos(t198));
	t225 = pkin(2) * qJD(2);
	t224 = t193 * t197;
	t199 = sin(pkin(8));
	t222 = t199 * t201;
	t221 = t199 * t202;
	t200 = cos(pkin(8));
	t220 = t200 * t201;
	t219 = t200 * t202;
	t216 = t227 * t199;
	t215 = t227 * t200;
	t208 = r_i_i_C(1) * t201 + r_i_i_C(2) * t202;
	t207 = t192 * t197 * t208;
	t205 = -cos(qJ(2)) * t225 + t204;
	t203 = t192 * r_i_i_C(2) * t223 - t208 * t193 * qJD(5) + (-pkin(3) * sin(t198) + t214 * t192) * t197 + t228 * t224;
	t1 = [0, t205 * t200 + t215, t200 * t204 + t215, 0, t200 * t207 + ((-t193 * t219 - t222) * r_i_i_C(1) + (t193 * t220 - t221) * r_i_i_C(2)) * qJD(5); 0, t205 * t199 + t216, t199 * t204 + t216, 0, t199 * t207 + ((-t193 * t221 + t220) * r_i_i_C(1) + (t193 * t222 + t219) * r_i_i_C(2)) * qJD(5); 0, -sin(qJ(2)) * t225 + t203, t203, 0, (-t202 * t224 + t213) * r_i_i_C(2) - t229 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end