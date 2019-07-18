% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S4PRRR2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S4PRRR2_jacobig_rot_0_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_jacobig_rot_0_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_jacobig_rot_0_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_0_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:29
% EndTime: 2019-07-18 13:27:29
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
Jg_rot  = t1;
