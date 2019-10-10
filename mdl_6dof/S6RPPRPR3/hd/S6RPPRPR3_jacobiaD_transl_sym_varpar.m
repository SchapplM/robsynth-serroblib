% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(9);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0, 0, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (53->19), mult. (76->33), div. (0->0), fcn. (48->6), ass. (0->14)
	t24 = sin(qJ(4));
	t25 = cos(qJ(4));
	t34 = (r_i_i_C(1) * t25 - r_i_i_C(2) * t24) * qJD(4);
	t33 = qJD(1) * t24;
	t32 = qJD(1) * t25;
	t31 = qJD(4) * t24;
	t30 = qJD(4) * t25;
	t29 = -pkin(2) - pkin(7) - r_i_i_C(3);
	t27 = r_i_i_C(1) * t24 + r_i_i_C(2) * t25 + qJ(3);
	t26 = qJD(3) + t34;
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t26 * t22 + (-cos(qJ(1)) * pkin(1) + t29 * t22 - t27 * t21) * qJD(1), 0, qJD(1) * t22, (-t21 * t30 - t22 * t33) * r_i_i_C(2) + (-t21 * t31 + t22 * t32) * r_i_i_C(1), 0, 0; t26 * t21 + (-sin(qJ(1)) * pkin(1) + t27 * t22 + t29 * t21) * qJD(1), 0, qJD(1) * t21, (-t21 * t33 + t22 * t30) * r_i_i_C(2) + (t21 * t32 + t22 * t31) * r_i_i_C(1), 0, 0; 0, 0, 0, -t34, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (97->21), mult. (106->27), div. (0->0), fcn. (69->8), ass. (0->17)
	t31 = qJ(4) + pkin(10);
	t27 = sin(t31);
	t29 = cos(t31);
	t40 = sin(qJ(4)) * pkin(4) + r_i_i_C(1) * t27 + r_i_i_C(2) * t29;
	t44 = qJD(4) * t40;
	t39 = cos(qJ(4)) * pkin(4) + r_i_i_C(1) * t29 - r_i_i_C(2) * t27;
	t43 = t39 * qJD(4);
	t32 = qJ(1) + pkin(9);
	t28 = sin(t32);
	t42 = qJD(1) * t28;
	t41 = -pkin(2) - r_i_i_C(3) - qJ(5) - pkin(7);
	t38 = qJ(3) + t40;
	t37 = qJD(1) * t39;
	t36 = qJD(3) + t43;
	t30 = cos(t32);
	t26 = qJD(1) * t30;
	t1 = [-t28 * qJD(5) + t36 * t30 + (-cos(qJ(1)) * pkin(1) + t41 * t30 - t38 * t28) * qJD(1), 0, t26, -t28 * t44 + t30 * t37, -t42, 0; t30 * qJD(5) + t36 * t28 + (-sin(qJ(1)) * pkin(1) + t41 * t28 + t38 * t30) * qJD(1), 0, t42, t28 * t37 + t30 * t44, t26, 0; 0, 0, 0, -t43, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:44
	% EndTime: 2019-10-09 23:37:44
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (272->46), mult. (312->71), div. (0->0), fcn. (240->10), ass. (0->36)
	t220 = qJ(4) + pkin(10);
	t216 = sin(t220);
	t218 = cos(t220);
	t247 = pkin(8) + r_i_i_C(3);
	t234 = t247 * t218 - sin(qJ(4)) * pkin(4);
	t223 = sin(qJ(6));
	t225 = cos(qJ(6));
	t235 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(5);
	t256 = (-t235 * t216 + t234) * qJD(4);
	t255 = -pkin(5) * t216 - qJ(3) + t234;
	t232 = t247 * t216 + cos(qJ(4)) * pkin(4);
	t253 = t235 * t218 + t232;
	t237 = qJD(1) * t216 + qJD(6);
	t252 = t223 * t237;
	t251 = t225 * t237;
	t243 = -pkin(2) - qJ(5) - pkin(7);
	t221 = qJ(1) + pkin(9);
	t217 = sin(t221);
	t242 = qJD(1) * t217;
	t241 = qJD(4) * t223;
	t240 = qJD(4) * t225;
	t239 = qJD(6) * t218;
	t238 = -qJD(6) * t216 - qJD(1);
	t236 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
	t231 = qJD(6) * t236;
	t230 = -t218 * t241 + t238 * t225;
	t229 = t218 * t240 + t238 * t223;
	t228 = qJD(3) + (pkin(5) * t218 + t232) * qJD(4);
	t227 = qJD(1) * t253;
	t219 = cos(t221);
	t215 = qJD(1) * t219;
	t214 = t229 * t217 + t219 * t251;
	t213 = t230 * t217 - t219 * t252;
	t212 = -t217 * t251 + t229 * t219;
	t211 = t217 * t252 + t230 * t219;
	t1 = [t212 * r_i_i_C(1) + t211 * r_i_i_C(2) - t217 * qJD(5) + t228 * t219 + (-cos(qJ(1)) * pkin(1) + t243 * t219 + t255 * t217) * qJD(1), 0, t215, t219 * t227 + (-t236 * t239 + t256) * t217, -t242, t213 * r_i_i_C(1) - t214 * r_i_i_C(2); t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t219 * qJD(5) + t228 * t217 + (-sin(qJ(1)) * pkin(1) + t243 * t217 - t255 * t219) * qJD(1), 0, t242, t217 * t227 + (t218 * t231 - t256) * t219, t215, -t211 * r_i_i_C(1) + t212 * r_i_i_C(2); 0, 0, 0, -t253 * qJD(4) + t216 * t231, 0, (t216 * t240 + t223 * t239) * r_i_i_C(2) + (t216 * t241 - t225 * t239) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end