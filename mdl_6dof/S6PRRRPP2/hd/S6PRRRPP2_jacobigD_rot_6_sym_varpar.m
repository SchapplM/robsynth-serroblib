% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPP2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:33
% EndTime: 2019-02-26 20:09:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t144 = sin(pkin(6));
t147 = sin(qJ(3));
t157 = t144 * t147;
t148 = sin(qJ(2));
t156 = t144 * t148;
t146 = cos(pkin(6));
t155 = t146 * t148;
t150 = cos(qJ(2));
t154 = t146 * t150;
t153 = qJD(2) * t147;
t143 = sin(pkin(10));
t145 = cos(pkin(10));
t152 = t143 * t150 + t145 * t155;
t151 = -t143 * t155 + t145 * t150;
t149 = cos(qJ(3));
t1 = [0, 0, t151 * qJD(2) (t143 * t157 + t151 * t149) * qJD(3) + (-t143 * t154 - t145 * t148) * t153, 0, 0; 0, 0, t152 * qJD(2) (-t145 * t157 + t152 * t149) * qJD(3) + (-t143 * t148 + t145 * t154) * t153, 0, 0; 0, 0, qJD(2) * t156, t144 * t150 * t153 + (t146 * t147 + t149 * t156) * qJD(3), 0, 0;];
JgD_rot  = t1;
