% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PRRRR1_jacobig_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobig_rot_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobig_rot_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.01s
% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
t1 = [0, 0, sin(qJ(2)), 0, 0; 0, -1, 0, 0, 0; 0, 0, -cos(qJ(2)), 0, 0;];
Jg_rot  = t1;
