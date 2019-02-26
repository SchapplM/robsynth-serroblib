% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRRP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRRP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRRP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:37:00
% EndTime: 2019-02-26 19:37:00
% DurationCPUTime: 0.05s
% Computational Cost: add. (98->11), mult. (48->12), div. (0->0), fcn. (24->6), ass. (0->14)
t58 = -pkin(3) - r_i_i_C(1);
t50 = qJD(1) + qJD(2);
t57 = pkin(2) * t50;
t56 = pkin(1) * qJD(1);
t51 = qJ(1) + qJ(2);
t49 = qJ(3) + t51;
t44 = sin(t49);
t45 = cos(t49);
t46 = qJD(3) + t50;
t55 = (r_i_i_C(2) * t44 + t45 * t58) * t46;
t54 = (-r_i_i_C(2) * t45 + t58 * t44) * t46;
t53 = -cos(t51) * t57 + t55;
t52 = -sin(t51) * t57 + t54;
t1 = [-cos(qJ(1)) * t56 + t53, t53, t55, 0; -sin(qJ(1)) * t56 + t52, t52, t54, 0; 0, 0, 0, 0;];
JaD_transl  = t1;
