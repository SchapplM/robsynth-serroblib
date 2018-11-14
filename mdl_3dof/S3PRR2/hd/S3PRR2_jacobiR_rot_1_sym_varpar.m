% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S3PRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
%
% Output:
% JR_rot [9x3]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S3PRR2_jacobiR_rot_1_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_jacobiR_rot_1_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_jacobiR_rot_1_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_1_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:57
% EndTime: 2018-11-14 10:12:57
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0;];
JR_rot  = t1;
