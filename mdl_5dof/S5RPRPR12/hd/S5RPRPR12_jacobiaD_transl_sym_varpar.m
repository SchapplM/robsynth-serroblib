% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:05
	% EndTime: 2019-12-31 18:31:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:05
	% EndTime: 2019-12-31 18:31:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(6) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(8) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(8)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t28, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:05
	% EndTime: 2019-12-31 18:31:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (116->28), mult. (184->43), div. (0->0), fcn. (139->7), ass. (0->18)
	t172 = pkin(8) + qJ(3);
	t170 = sin(t172);
	t171 = cos(t172);
	t173 = sin(pkin(9));
	t174 = cos(pkin(9));
	t182 = r_i_i_C(1) * t174 - r_i_i_C(2) * t173 + pkin(3);
	t188 = r_i_i_C(3) + qJ(4);
	t190 = (t182 * t170 - t188 * t171) * qJD(3) - t170 * qJD(4);
	t176 = sin(qJ(1));
	t187 = qJD(1) * t176;
	t177 = cos(qJ(1));
	t186 = qJD(1) * t177;
	t185 = qJD(3) * t171;
	t183 = qJD(3) * t188;
	t181 = t173 * r_i_i_C(1) + t174 * r_i_i_C(2) + pkin(6) + qJ(2);
	t180 = -t182 * qJD(3) + qJD(4);
	t179 = -t188 * t170 - t182 * t171 - cos(pkin(8)) * pkin(2) - pkin(1);
	t1 = [t177 * qJD(2) + t190 * t176 + (-t181 * t176 + t179 * t177) * qJD(1), t186, (-t177 * t183 + t182 * t187) * t170 + (t180 * t177 - t188 * t187) * t171, -t170 * t187 + t177 * t185, 0; t176 * qJD(2) - t190 * t177 + (t179 * t176 + t181 * t177) * qJD(1), t187, (-t176 * t183 - t182 * t186) * t170 + (t180 * t176 + t188 * t186) * t171, t170 * t186 + t176 * t185, 0; 0, 0, -t190, qJD(3) * t170, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:05
	% EndTime: 2019-12-31 18:31:06
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (257->47), mult. (309->75), div. (0->0), fcn. (248->9), ass. (0->38)
	t212 = cos(pkin(9)) * pkin(4) + pkin(3);
	t219 = pkin(8) + qJ(3);
	t215 = sin(t219);
	t239 = t215 * qJD(4);
	t217 = cos(t219);
	t248 = r_i_i_C(3) + pkin(7) + qJ(4);
	t250 = t248 * t217;
	t254 = (-t212 * t215 + t250) * qJD(3) + t239;
	t218 = pkin(9) + qJ(5);
	t214 = sin(t218);
	t216 = cos(t218);
	t229 = r_i_i_C(1) * t216 - r_i_i_C(2) * t214 + t212;
	t226 = -t229 * t215 + t250;
	t224 = cos(qJ(1));
	t240 = qJD(5) * t217;
	t233 = -qJD(1) + t240;
	t252 = t224 * t233;
	t232 = qJD(1) * t217 - qJD(5);
	t223 = sin(qJ(1));
	t243 = qJD(3) * t223;
	t249 = -t215 * t243 + t232 * t224;
	t246 = qJD(1) * t223;
	t245 = qJD(1) * t224;
	t244 = qJD(3) * t217;
	t242 = qJD(3) * t224;
	t241 = qJD(5) * t215;
	t237 = t248 * t215;
	t234 = pkin(4) * sin(pkin(9)) + pkin(6) + qJ(2);
	t231 = r_i_i_C(1) * t214 + r_i_i_C(2) * t216;
	t230 = t233 * t223;
	t228 = -t212 * t217 - cos(pkin(8)) * pkin(2) - pkin(1) - t237;
	t227 = t215 * t242 + t232 * t223;
	t225 = qJD(4) * t217 + t231 * t241 + (-t229 * t217 - t237) * qJD(3);
	t211 = t214 * t230 - t249 * t216;
	t210 = t249 * t214 + t216 * t230;
	t209 = t214 * t252 + t227 * t216;
	t208 = t227 * t214 - t216 * t252;
	t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t224 * qJD(2) - t254 * t223 + (-t234 * t223 + t228 * t224) * qJD(1), t245, t225 * t224 - t226 * t246, -t215 * t246 + t217 * t242, t208 * r_i_i_C(1) + t209 * r_i_i_C(2); -t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t223 * qJD(2) + t254 * t224 + (t228 * t223 + t234 * t224) * qJD(1), t246, t225 * t223 + t226 * t245, t215 * t245 + t217 * t243, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2); 0, 0, t226 * qJD(3) - t231 * t240 + t239, qJD(3) * t215, (t214 * t241 - t216 * t244) * r_i_i_C(2) + (-t214 * t244 - t216 * t241) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end