% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:22
% EndTime: 2019-02-26 21:25:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (103->20), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->17)
t149 = qJ(2) + pkin(9);
t147 = sin(t149);
t148 = cos(t149);
t165 = r_i_i_C(3) + qJ(4);
t169 = pkin(3) - r_i_i_C(2);
t158 = t169 * t147 - t165 * t148 + sin(qJ(2)) * pkin(2);
t155 = -t158 * qJD(2) + t147 * qJD(4);
t172 = t155 + (r_i_i_C(1) + qJ(3) + pkin(7)) * qJD(1);
t171 = -t165 * t147 - t169 * t148 - cos(qJ(2)) * pkin(2);
t152 = sin(qJ(1));
t164 = qJD(1) * t152;
t154 = cos(qJ(1));
t163 = qJD(1) * t154;
t162 = qJD(2) * t148;
t157 = qJD(3) + (-pkin(1) + t171) * qJD(1);
t156 = t171 * qJD(2) + qJD(4) * t148;
t1 = [-t172 * t152 + t157 * t154, t156 * t154 + t158 * t164, t163, -t147 * t164 + t154 * t162, 0, 0; t157 * t152 + t172 * t154, t156 * t152 - t158 * t163, t164, t147 * t163 + t152 * t162, 0, 0; 0, t155, 0, qJD(2) * t147, 0, 0;];
JaD_transl  = t1;
