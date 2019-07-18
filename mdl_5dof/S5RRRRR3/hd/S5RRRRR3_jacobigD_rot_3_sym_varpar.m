% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRRR3_jacobigD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobigD_rot_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobigD_rot_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobigD_rot_3_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
t27 = qJD(1) * cos(qJ(1));
t26 = qJD(1) * sin(qJ(1));
t1 = [0, t27, t27, 0, 0; 0, t26, t26, 0, 0; 0, 0, 0, 0, 0;];
JgD_rot  = t1;
