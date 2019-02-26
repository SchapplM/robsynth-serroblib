% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR9_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobigD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:28
% EndTime: 2019-02-26 22:20:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
t148 = sin(pkin(6));
t161 = t148 * cos(pkin(7));
t151 = sin(qJ(2));
t152 = sin(qJ(1));
t160 = t151 * t152;
t154 = cos(qJ(1));
t159 = t151 * t154;
t153 = cos(qJ(2));
t158 = t152 * t153;
t157 = t153 * t154;
t156 = qJD(1) * t148;
t147 = sin(pkin(7));
t155 = qJD(2) * t147;
t150 = cos(pkin(6));
t1 = [0, t154 * t156 -(t150 * t160 - t157) * t155 + (-(-t150 * t157 + t160) * t147 + t154 * t161) * qJD(1), 0, 0, 0; 0, t152 * t156 -(-t150 * t159 - t158) * t155 + (-(-t150 * t158 - t159) * t147 + t152 * t161) * qJD(1), 0, 0, 0; 0, 0, t148 * t151 * t155, 0, 0, 0;];
JgD_rot  = t1;
