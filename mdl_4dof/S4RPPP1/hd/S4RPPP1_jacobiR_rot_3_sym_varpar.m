% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S4RPPP1_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiR_rot_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiR_rot_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->8), mult. (18->14), div. (0->0), fcn. (24->9), ass. (0->10)
t56 = cos(qJ(1));
t55 = sin(qJ(1));
t54 = cos(pkin(6));
t53 = sin(pkin(4));
t52 = sin(pkin(6));
t51 = pkin(4) - pkin(6);
t50 = pkin(4) + pkin(6);
t49 = cos(t51) / 0.2e1 + cos(t50) / 0.2e1;
t48 = sin(t50) / 0.2e1 - sin(t51) / 0.2e1;
t1 = [t56 * t53, 0, 0, 0; t55 * t53, 0, 0, 0; 0, 0, 0, 0; t56 * t48 + t55 * t54, 0, 0, 0; t55 * t48 - t56 * t54, 0, 0, 0; 0, 0, 0, 0; t56 * t49 - t55 * t52, 0, 0, 0; t55 * t49 + t56 * t52, 0, 0, 0; 0, 0, 0, 0;];
JR_rot  = t1;
